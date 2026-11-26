# LaTeX, Markdown & Biomedical Typography Reference

Load this reference ONLY when cleaning up LaTeX/Markdown layout, journal typography, or symbol formatting.

## 1. Syntax & Anchor Preservation Rules

- **Never break**: `\cite{}`, `\citep{}`, `\citet{}`, `\ref{}`, `\eqref{}`, `\label{}`, BibTeX keys, Pandoc citations (`[@key]`), equation environments, or Markdown table pipe counts.
- **Non-breaking ties**: Keep `Fig.~\ref{...}`, `Table~\ref{...}`, `Extended Data Fig.~\ref{...}`, and `5~mg/mL` in LaTeX to prevent orphan line breaks.
- **Dashes**: Use hyphen `-` for compound modifiers (`single-cell`), en-dash `--` (`–`) for ranges (`10--20`, `Fig. 2A--C`), and em-dash `---` (`—`) sparingly.

## 2. Biomedical & Statistical Typography

| Item | Convention | Example |
| --- | --- | --- |
| Statistical symbols | Italicized in formatted prose/LaTeX | *P* = 0.004, *n* = 32, *r* = 0.61, *t*, *F* |
| Human gene vs protein | Gene: ALL CAPS *italic*; Protein: ALL CAPS roman | *TP53* gene; TP53 protein |
| Mouse gene vs protein | Gene: Title Case *italic*; Protein: ALL CAPS roman | *Trp53* gene; TRP53 protein |
| Species & Latin terms | *Italicized* | *Homo sapiens*, *Mus musculus*, *in vitro*, *in vivo*, *de novo* |
| Units & spaces | Space between value and unit (except `%` and `°C` per journal) | `5 mg/mL`, `37 °C`, `24.3%`, `P < 0.001` |

## 3. Tense Consistency Quick Check

- **Past tense**: Specific Methods you executed and Results you observed (*"We sequenced..."*, *"Expression increased..."*).
- **Present tense**: Established literature facts, mathematical definitions, and direct pointers to figures/tables (*"Fig. 2A shows..."*, *"Cox regression models..."*).
