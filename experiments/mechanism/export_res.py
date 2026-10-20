"""Pack light copies of /data/results_mechanism (json + R / light arrays) into out/export.tar.gz."""
import os, glob, tarfile, io, json
import numpy as np
RES = os.environ.get("RES_DIR", "/data/results_mechanism")  # Modal volume path; override locally; os.makedirs("out", exist_ok=True)
n = 0
with tarfile.open("out/mech_r3_export.tar.gz", "w:gz") as t:
    for f in glob.glob(RES + "/*.npz"):
        t.add(f, arcname="res/" + os.path.basename(f))
    for arm_dir in sorted(glob.glob(RES + "/*")):
        if not os.path.isdir(arm_dir):
            continue
        arm = os.path.basename(arm_dir)
        for j in glob.glob(arm_dir + "/*.json"):
            cid = os.path.basename(j)[:-5]; t.add(j, arcname=f"res/{arm}/{cid}.json"); n += 1
            for ext in (".anat.npz", ".npz"):
                p = os.path.join(arm_dir, cid + ext)
                if not os.path.exists(p):
                    continue
                if arm in ("SN", "N", "SY", "Y") or ext == ".anat.npz":
                    t.add(p, arcname=f"res/{arm}/{cid}{ext}")
                else:
                    z = np.load(p); keep = {k: z[k] for k in z.files if k == "R" or k.startswith("trace_")}
                    b = io.BytesIO(); np.savez_compressed(b, **keep); b.seek(0)
                    ti = tarfile.TarInfo(f"res/{arm}/{cid}.npz"); ti.size = len(b.getvalue()); t.addfile(ti, b)
                break
print("exported", n, "conditions")
