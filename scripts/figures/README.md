# Supplementary Figures With Caption

All supplementary figures + captions in the same PDF file.

Tips:
- avoid dots in the figure names (LaTeX doen's like them sometimes)

Useful instructions:

```shell
pdfjam --nup 3x8 SupplementaryFigure36.pggb.wgg.88.bisnips.chm13_vs_grch38.pdf 
pdfcrop --margins "1 1 1 1" SupplementaryFigure36.pggb.wgg.88.bisnips.chm13_vs_grch38-pdfjam.pdf SupplementaryFigure36.pdf
```
