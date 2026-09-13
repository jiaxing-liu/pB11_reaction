# Independent fusion library program

Current approved work is coordinated by /home/cloud/research/pB-baldur/baldur-code/docs/exl50u-program/PLAN.md and STATUS.md. Read these and the relevant API_CONTRACT.md before continuing; check actual git state after resuming.

Work on feature/exl50u-fusion-library. Preserve the existing scalar C++/C/Fortran interfaces; use explicit units and model/domain contracts. Keep BALDUR-specific state, indexing and I/O out of this library. Root agent owns physics and architecture; Luna may implement bounded agreed interfaces/tests. Never silently add empirical physics or turn failed calculations into zero-valued success.

Commit and push reviewed, validated milestones to this same branch; do not merge main. Record paired host/library commits in the coordinator STATUS.md. Do not overwrite concurrent agents' owned files.
