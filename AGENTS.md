# Independent fusion library program

Current approved work is coordinated by /home/cloud/research/pB-baldur/baldur-code/docs/exl50u-program/PLAN.md and STATUS.md. Read these and the relevant API_CONTRACT.md before continuing; check actual git state after resuming.

Work on feature/exl50u-fusion-library. Preserve the existing scalar C++/C/Fortran interfaces; use explicit units and model/domain contracts. Keep BALDUR-specific state, indexing and I/O out of this library. Root agent owns physics and architecture; Luna may implement bounded agreed interfaces/tests. Never silently add empirical physics or turn failed calculations into zero-valued success.

Commit and push reviewed, validated milestones to this same branch; do not merge main. Record paired host/library commits in the coordinator STATUS.md. Do not overwrite concurrent agents' owned files.

## Chat-first analysis routing — confirmed 2026-09-20

Follow the coordinating BALDUR AGENTS.md cost preference across sessions.
Route physics reasoning, design review, diagnosis, test planning and result
interpretation to ordinary Chat whenever practical. Codex collects focused
project evidence, implements and verifies; avoid duplicate extended analysis
and new Codex analysis agents. The authorized preferred Chat is `讨论pB聚变`
(id `6aae7d5a-22f0-83ea-bd79-de8842487e9c`), called “讨论pB反应” by the user.
Send only necessary project materials; record adopted advice in project docs.
Do not silently substitute Work/API/Codex analysis or claim free usage/Pro
selection without verification. Continue independent execution while awaiting
Chat; record unavailability if it occurs.
