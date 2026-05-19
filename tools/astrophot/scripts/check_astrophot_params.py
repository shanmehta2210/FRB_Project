"""Quick check of AstroPhot parameter names, locking API, chi2 computation."""
import astrophot as ap
import numpy as np

# Create a dummy model to see parameter names
target = ap.image.Target_Image(data=np.random.randn(50,50), pixelscale=0.258)
model = ap.models.AstroPhot_Model(
    name='test',
    model_type='sersic galaxy model',
    target=target,
)
model.initialize()

# List all parameters
print("=== Parameters ===")
for p in model.parameter_order:
    param = model[p]
    print(f"  {p}: value={param.value}, locked={param.locked}")

print("\nParameter names:", model.parameter_order)

# Check locking API
print("\n=== Locking test ===")
model["q"].value = 0.5
model["q"].locked = True
print(f"q locked={model['q'].locked}, value={model['q'].value}")

# Unlock q for fitting test
model["q"].locked = False

# Run a fit first so we have loss history
fitter = ap.fit.LM(model, verbose=0)
result = fitter.fit()
loss = fitter.res_loss()
ndf = fitter.ndf
print(f"\nres_loss = {loss}")
print(f"ndf = {ndf}")
print(f"loss/ndf = {loss/ndf}")
print(f"Type of loss: {type(loss)}")

# Compute chi2 manually from model image
model_image = model()
residuals = (target.data.detach().numpy() - model_image.data.detach().numpy())
var = target.variance.detach().numpy()
chi2_manual = np.sum(residuals**2 / var)
n_pixels = residuals.size
n_params = len([p for p in model.parameter_order if not model[p].locked])
ndf_manual = n_pixels - n_params
chi2nu_manual = chi2_manual / ndf_manual
print(f"\nManual chi2_total = {chi2_manual}")
print(f"Manual ndf = {ndf_manual}")
print(f"Manual chi2_nu = {chi2nu_manual}")
print(f"\nComparison:")
print(f"  res_loss vs chi2_nu:  {loss:.6f} vs {chi2nu_manual:.6f}")
print(f"  res_loss vs chi2:     {loss:.6f} vs {chi2_manual:.6f}")
print(f"  res_loss * ndf =      {loss * ndf:.6f}")
print(f"  Ratio res_loss/chi2_nu = {loss/chi2nu_manual:.6f}")
print(f"  Ratio res_loss/chi2    = {loss/chi2_manual:.6f}")

# Also check: is loss_history storing chi2 or chi2_nu?
print(f"\nLoss history (last 3): {fitter.loss_history[-3:]}")

# Try getting model's own loss
import torch
model_img = model()
diff = target.data - model_img.data
weighted_resid = diff / torch.sqrt(target.variance)
chi2_torch = torch.sum(weighted_resid**2).item()
print(f"\nTorch manual chi2 = {chi2_torch}")
print(f"Torch chi2/ndf = {chi2_torch / ndf}")
