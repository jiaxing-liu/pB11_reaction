# Independent fusion library program

Current approved work is coordinated by /home/cloud/research/pB-baldur/baldur-code/docs/exl50u-program/PLAN.md and STATUS.md. Read these and the relevant API_CONTRACT.md before continuing; check actual git state after resuming.

Work on feature/exl50u-fusion-library. Preserve the existing scalar C++/C/Fortran interfaces; use explicit units and model/domain contracts. Keep BALDUR-specific state, indexing and I/O out of this library. Root agent owns physics and architecture; bounded agreed tasks may be delegated to gpt-6.1-sol. Never silently add empirical physics or turn failed calculations into zero-valued success.

Commit and push reviewed, validated milestones to this same branch; do not merge main. Record paired host/library commits in the coordinator STATUS.md. Do not overwrite concurrent agents' owned files.

## Current model rule — confirmed 2026-10-07

The latest user instruction requires `gpt-6.1-sol` for every task, including
analysis, implementation, execution, documentation and review. Explicitly select
`gpt-6.1-sol` for all new subagents. Do not assign or resume work on existing
`gpt-6-sol` or `gpt-6-luna` agents. The user selects the root model in the client;
tools cannot change it. This rule supersedes all historical Chat/Pro/Luna/sol
routing below and in older program records. Follow the current coordinating
BALDUR AGENTS.md and preserve its case-concurrency and memory limits.

## Historical Chat-first routing — superseded 2026-10-07

The following destination and authorization are retained as historical records;
they are not the current task-routing policy.

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
