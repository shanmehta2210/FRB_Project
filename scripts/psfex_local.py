import numpy as np
from astropy.io import fits


class PSFEx:
    """Lightweight reader/evaluator for PSFEx .psf models.

    Provides the minimal interface used in this project:
    `PSFEx(path).get_rec(y, x)`.
    """

    def __init__(self, psf_path):
        with fits.open(psf_path) as hdul:
            hdr = hdul[1].header
            row = hdul[1].data[0]

        self.polnaxis = int(hdr.get("POLNAXIS", 2))
        self.poldeg = int(hdr.get("POLDEG1", 0))
        self.psf_samp = float(hdr.get("PSF_SAMP", 1.0))
        if self.psf_samp <= 0.0:
            raise ValueError(f"Invalid PSF_SAMP={self.psf_samp}")
        self.pixstep = 1.0 / self.psf_samp

        self.polzero = [float(hdr.get(f"POLZERO{i+1}", 0.0)) for i in range(self.polnaxis)]
        self.polscal = [float(hdr.get(f"POLSCAL{i+1}", 1.0)) for i in range(self.polnaxis)]

        # Native PSF basis cube: [ncoeff, ny, nx]
        self.psf_mask = np.asarray(row["PSF_MASK"], dtype=np.float64)

        if self.psf_mask.ndim != 3:
            raise ValueError(f"Unexpected PSF_MASK shape: {self.psf_mask.shape}")

        self.ncoeff = self.psf_mask.shape[0]
        self.mask_ny = int(self.psf_mask.shape[1])
        self.mask_nx = int(self.psf_mask.shape[2])

        # Upstream wrapper expects 2 context dimensions.
        if self.polnaxis != 2:
            raise ValueError(f"Expected POLNAXIS=2, got {self.polnaxis}")

        # Match upstream reconstruction-size rule: ceil(mask_dim*psf_samp), forced odd.
        self.recon_nx = int(np.ceil(self.mask_nx * self.psf_samp))
        self.recon_ny = int(np.ceil(self.mask_ny * self.psf_samp))
        if (self.recon_nx % 2) == 0:
            self.recon_nx += 1
        if (self.recon_ny % 2) == 0:
            self.recon_ny += 1

        self.exponents = self._build_exponents_upstream_2d(self.poldeg)

        if len(self.exponents) != self.ncoeff:
            raise ValueError(
                f"Coefficient count mismatch: model has {self.ncoeff}, expected {len(self.exponents)} "
                f"for POLNAXIS={self.polnaxis}, POLDEG={self.poldeg}"
            )

    @staticmethod
    def _build_exponents_upstream_2d(deg):
        # Mirrors poly.c term traversal for POLY_DIM=2, single group.
        # Order: (0,0), (1,0), (2,0), ..., then terms with y-power 1, etc.
        exps = [(0, 0)]
        for ypow in range(0, deg + 1):
            x_start = 1 if ypow == 0 else 0
            for xpow in range(x_start, deg - ypow + 1):
                exps.append((xpow, ypow))
        return exps

    @staticmethod
    def _interpf(x):
        # Lanczos kernel used in upstream psfex.h (INTERPF, INTERPFAC=4).
        if -1.0e-5 < x < 1.0e-5:
            return 1.0
        if x > 4.0 or x < -4.0:
            return 0.0
        pix = np.pi * x
        return np.sin(pix) * np.sin(0.25 * pix) / (0.25 * (pix * pix))

    def _vignet_resample(self, pix1, w2, h2, dx, dy, step2):
        # Python port of upstream _psfex_vignet_resample from src/psfex.c.
        h1, w1 = pix1.shape
        pix2 = np.zeros((h2, w2), dtype=np.float64)

        xc1 = float(w1 // 2)
        xc2 = float(w2 // 2)
        xs1 = xc1 + dx - xc2 * step2

        if int(xs1) >= w1:
            return pix2

        ixs2 = 0
        if xs1 < 0.0:
            dix2 = int(1 - xs1 / step2)
            if dix2 >= w2:
                return pix2
            ixs2 += dix2
            xs1 += dix2 * step2

        nx2 = int((w1 - 1 - xs1) / step2 + 1)
        ix2 = w2 - ixs2
        if nx2 > ix2:
            nx2 = ix2
        if nx2 <= 0:
            return pix2

        yc1 = float(h1 // 2)
        yc2 = float(h2 // 2)
        ys1 = yc1 + dy - yc2 * step2
        if int(ys1) >= h1:
            return pix2

        iys2 = 0
        if ys1 < 0.0:
            diy2 = int(1 - ys1 / step2)
            if diy2 >= h2:
                return pix2
            iys2 += diy2
            ys1 += diy2 * step2

        ny2 = int((h1 - 1 - ys1) / step2 + 1)
        iy2 = h2 - iys2
        if ny2 > iy2:
            ny2 = iy2
        if ny2 <= 0:
            return pix2

        iys1a = int(ys1)
        hmh = (8 // 2) - 1
        if iys1a < 0 or ((iys1a - hmh) < 0):
            iys1a = 0
        else:
            iys1a = iys1a - hmh

        ny1 = int(ys1 + ny2 * step2) + 8 - hmh
        if ny1 > h1:
            ny1 = h1
        ny1 -= iys1a
        ys1 -= float(iys1a)

        # X-direction local interpolation masks.
        hmw = (8 // 2) - 1
        x_starts = np.zeros(nx2, dtype=np.int64)
        x_weights = []
        x1 = xs1
        for j in range(nx2):
            ix1 = int(x1)
            ix = ix1 - hmw
            dxm = ix1 - x1 - hmw
            if ix < 0:
                n = 8 + ix
                dxm -= float(ix)
                ix = 0
            else:
                n = 8
            t = w1 - ix
            if n > t:
                n = t

            x_starts[j] = ix
            w = np.empty(n, dtype=np.float64)
            xx = dxm
            for i in range(n):
                w[i] = self._interpf(xx)
                xx += 1.0
            s = float(np.sum(w))
            if s > 0.0:
                w /= s
            x_weights.append(w)
            x1 += step2

        # Interpolate in x and transpose into [nx2, ny1].
        pix12 = np.zeros((nx2, ny1), dtype=np.float64)
        for k in range(ny1):
            row_src = pix1[iys1a + k, :]
            for j in range(nx2):
                st = int(x_starts[j])
                w = x_weights[j]
                pix12[j, k] = np.dot(w, row_src[st : st + w.size])

        # Y-direction local interpolation masks.
        y_starts = np.zeros(ny2, dtype=np.int64)
        y_weights = []
        y1 = ys1
        for j in range(ny2):
            iy1 = int(y1)
            iy = iy1 - hmh
            dym = iy1 - y1 - hmh
            if iy < 0:
                n = 8 + iy
                dym -= float(iy)
                iy = 0
            else:
                n = 8
            t = ny1 - iy
            if n > t:
                n = t

            y_starts[j] = iy
            w = np.empty(n, dtype=np.float64)
            yy = dym
            for i in range(n):
                w[i] = self._interpf(yy)
                yy += 1.0
            s = float(np.sum(w))
            if s > 0.0:
                w /= s
            y_weights.append(w)
            y1 += step2

        # Interpolate in y and transpose back into output image.
        for k in range(nx2):
            col_src = pix12[k, :]
            out_x = ixs2 + k
            for j in range(ny2):
                st = int(y_starts[j])
                w = y_weights[j]
                out_y = iys2 + j
                pix2[out_y, out_x] = np.dot(w, col_src[st : st + w.size])

        return pix2

    def get_rec(self, y, x, out_size=None):
        coords = [float(x), float(y)]

        # Normalize coordinates with PSFEx polynomial scaling.
        u = []
        for i in range(self.polnaxis):
            scale = self.polscal[i] if self.polscal[i] != 0 else 1.0
            u.append((coords[i] - self.polzero[i]) / scale)

        weights = np.empty(self.ncoeff, dtype=np.float64)
        for i, exp in enumerate(self.exponents):
            w = 1.0
            for axis, power in enumerate(exp):
                if power:
                    w *= u[axis] ** power
            weights[i] = w

        # Weighted combination in model grid.
        maskloc = np.tensordot(weights, self.psf_mask, axes=(0, 0))

        # Match upstream subpixel recentering and sampling step.
        dcol = float(x) - float(int(float(x) + 0.5))
        drow = float(y) - float(int(float(y) + 0.5))

        if out_size is None:
            out_nx = self.recon_nx
            out_ny = self.recon_ny
        else:
            out_nx = int(out_size)
            out_ny = int(out_size)
            if out_nx < 1 or out_ny < 1:
                raise ValueError(f"Invalid out_size={out_size}")
            if (out_nx % 2) == 0 or (out_ny % 2) == 0:
                raise ValueError(f"out_size must be odd, got {out_size}")

        rec = self._vignet_resample(
            maskloc,
            w2=out_nx,
            h2=out_ny,
            dx=-dcol * self.pixstep,
            dy=-drow * self.pixstep,
            step2=self.pixstep,
        )
        return rec.astype(np.float64)

    def get_center(self, y, x, out_size=None):
        if out_size is None:
            ncol = self.recon_nx
            nrow = self.recon_ny
        else:
            ncol = int(out_size)
            nrow = int(out_size)
        dcol = float(x) - float(int(float(x) + 0.5))
        drow = float(y) - float(int(float(y) + 0.5))
        colcen = float(ncol // 2) + dcol
        rowcen = float(nrow // 2) + drow
        return np.array([rowcen, colcen], dtype=np.float64)
