# Pushing to GitHub (not yet done)

The repository is committed locally on branch `main`. Nothing has been pushed. To publish it as a **private** repository under the KalarisLabs organisation:

```bash
cd fly-brain-resuscitation
# Option A: GitHub CLI (needs `gh auth login` with org repo-create rights)
gh repo create KalarisLabs/fly-brain-resuscitation --private --source=. --remote=origin \
  --description "In-silico anoxia, brain death and resuscitation on Drosophila connectomes (model predictions only)"
git push -u origin main

# Option B: create an empty private repo at https://github.com/organizations/KalarisLabs/repositories/new
#           (no README/licence/gitignore), then:
git remote add origin https://github.com/KalarisLabs/fly-brain-resuscitation.git
git push -u origin main
```

Before going public: confirm the FlyWire/Codex data terms in `DATA_LICENSE.md`, and check that `git ls-files` contains no raw data (`.gitignore` excludes `data/raw/`, `outputs/`, `*.npz` except the two driver rankings, `*.parquet`, archives). To make the repo public later: `gh repo edit KalarisLabs/fly-brain-resuscitation --visibility public --accept-visibility-change-consequences`.
