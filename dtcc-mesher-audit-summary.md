# DTCC Mesher vs Triangle — Copyright Audit Summary

Comparison of three independent audits examining whether `dtcc-mesher` infringes on Jonathan Shewchuk's Triangle library (https://www.cs.cmu.edu/~quake/triangle.html).

| Dimension | Claude (5 parallel agents) | Codex run #1 (7m 5s) | Codex run #2 (5m 36s) |
|---|---|---|---|
| **Verdict** | Clean | Clean | Clean |
| **Confidence** | High | High (explicit) | High (explicit) |
| **Scope** | 5 parallel passes: predicates, meshing fingerprints, global text, git history, public API | Single deep pass: predicates, meshing, style, residual-token grep, fixtures, build | Single deep pass: predicates, structures, algorithms, refinement control flow, API |
| **predicates.c check** | Public-domain header confirmed, 4262 lines of pure predicates | Same + noted 3D `orient3d`/`insphere` presence (stronger evidence of standalone release) | Same + explicitly cited CMU pages stating Triangle is copyrighted while robust predicates are public domain |
| **Triangle struct fingerprints (`otri`/`osub`/`triangulateio`/`mesh`/`behavior`)** | None found | None found | None found |
| **Triangle macro fingerprints (`bond`/`tspivot`/`sym`/`lnext`/`infect`)** | None found | None found | None found |
| **Triangle function fingerprints (`divconqdelaunay`/`plague`/`poolinit`)** | None found | None found | None found |
| **Architecture assessment** | Index-based (`TMTriangle{v[3], nbr[3]}`) vs Triangle's pointer-chasing — fundamentally different | Same, plus pairwise neighbor rebuild (`mesh.c:1158`) explicitly unlike Triangle's oriented-handle internals | Cited Shewchuk's paper §2 describing pointer-rich records — confirmed "not a rename-level rewrite" |
| **Algorithm selection** | Noted independent reimplementation | Noted incremental Bowyer-Watson with supertriangle | **Key insight: Triangle defaults to divide-and-conquer (`-i` is opt-in); dtcc-mesher uses incremental as only path → anti-fingerprint** |
| **Refinement scheduler** | Not analyzed at this depth | Candidate arrays + circumcenter insertion noted | **Key insight: Triangle uses 64 FIFO queues bucketed by angle; dtcc-mesher uses simple growable array with last-in pop → different engineering** |
| **Segment recovery** | Not analyzed at this depth | Flagged "especially non-Triangle-like": scan-crossings → flip-until-stall → rebuild-unconstrained fallback | Not flagged as deeply |
| **Acute-corner protection** | Not analyzed | Flagged custom `simple`/`shell` modes as distinctive to dtcc-mesher | Not flagged |
| **Git history** | 16 commits, single author, Mar 31–Apr 6 2026, no deleted Triangle files, pickaxe `-S` clean | Confirmed pickaxe results clean | Confirmed pickaxe results clean |
| **Public API** | Opaque `dtcc_mesher_domain`/`_mesh`; `num_*` naming; typed point/segment structs; long-form CLI flags | Flat index arrays, not `triangulateio` | Custom domain/mesh structs, `.pts`/`.pslg`/`.tri` workflow not Triangle's `.node`/`.ele`/`.poly` |
| **Test fixtures check** | Not checked for Triangle samples | **Explicitly verified absence of Triangle sample names** (`A.poly`, `box.poly`, `dots.poly`, `double.poly`, `face.poly`, `spiral.poly`) | Not checked |
| **Style/personality check** | Not analyzed | Noted absence of Shewchuk-style long prose / humor in core files; dtcc uses terse `tm_` prefix style | Noted core files start as "ordinary local C modules", not Shewchuk-style banner comments |
| **AI regurgitation risk** | Flagged as residual non-zero risk (README credits Codex) | Same; found no evidence of regurgitated islands | Same; noted cannot prove from repo alone that no one ever *looked* at `triangle.c` |
| **Novel false positives flagged and dismissed** | — | CMake description "Two-Dimensional Quality Mesh Generator" (`:3`) and log string "live subsegments" (`cdt.c:4456`) — both judged as generic terminology, not evidence | — |
| **Sources cited** | Internal forensic only | Internal forensic | Internal + 5 external CMU URLs (Triangle paper, help page, robust page) |
| **Recommendation** | Add upstream URL/checksum to `THIRD_PARTY.md` | Add checksum + note that `predicates.c` is the standalone release | Add retrieval date + checksum for commercial diligence |

## Where the three audits add distinct value

- **Claude (breadth)**: git history and global text scan — confirmed no Triangle code was ever present, even in deleted commits
- **Codex #1 (depth on style & fixtures)**: style-personality check, test fixture filename scan, segment-recovery + acute-protection specifics
- **Codex #2 (external comparison)**: grounded findings in Shewchuk's own published paper and docs — caught two *positive anti-fingerprints* (incremental-only triangulation, simple candidate array vs 64-queue scheduler) that materially strengthen the "independent implementation" conclusion

**Unanimous verdict: clean.** Three independent methodologies, zero Triangle meshing fingerprints, two positive anti-fingerprints from Codex #2, consistent public-domain attribution for the only vendored file.
