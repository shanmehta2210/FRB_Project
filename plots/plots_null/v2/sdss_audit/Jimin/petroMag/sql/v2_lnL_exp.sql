SELECT 
p.objid, p.ra, p.dec,
p.r AS modelMag_r,
p.cModelMag_u, p.cModelMag_g, p.cModelMag_r,
p.petroR90_r, p.petroR50_r,
p.petroMag_r, p.petroRad_r,
p.deVRad_r, p.deVAB_r, p.lnLDeV_r,
p.expRad_r, p.expAB_r, p.lnLExp_r,
p.lnLStar_r,
p.fracDeV_r, pz.nnAvgZ AS photz

FROM PhotoObj AS p
JOIN Photoz AS pz ON pz.objid = p.objid
WHERE
p.ra BETWEEN 148.0 AND 152.0
AND p.dec BETWEEN 0.0 AND 4.0
AND p.petroMag_r BETWEEN 12 AND 21
AND p.mode = 1
AND p.clean = 1
AND p.type_r = 3
AND p.lnLStar_r < -10
AND pz.nnAvgZ > 0
AND p.score > 0.8
AND p.lnLDeV_r < p.lnLExp_r
