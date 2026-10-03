# Colour-magnitude diagram of stars within 50 pc

`cmd_50pc.py` plots Tycho-2 B_T − V_T against absolute 2MASS K_s magnitude
(M_Ks = Ks + 5 log10(ϖ/mas) − 10) for Gaia DR3 stars with parallax > 20 mas.
It uses the Gaia DR3 best-neighbour cross-match tables.

```
pip install -r requirements.txt
python cmd_50pc.py            # writes cmd_50pc.png and caches data/*.csv
```

The script needs network access to `gea.esac.esa.int` (Gaia archive).
Quality cuts: ϖ/σϖ > 10, RUWE < 1.4, σ(B_T−V_T) < 0.1, 2MASS Ks quality "A".
Stars that fail the cuts are shown in grey.
