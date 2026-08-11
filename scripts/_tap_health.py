import time
import pyvo

t0 = time.time()
try:
    s = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")
    r = s.search("SELECT TOP 1 objid FROM ls_dr10.tractor")
    print("tap ok", len(r), f"{time.time()-t0:.1f}s")
except Exception as e:
    print("tap down", type(e).__name__, f"{time.time()-t0:.1f}s", str(e)[:180])
