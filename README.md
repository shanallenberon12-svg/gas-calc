# GasCalc Pro — Gas Compressibility Calculator Suite

A petroleum engineering web application with 5 calculators.

## Run Locally

```bash
pip install flask numpy
python app.py
```
Then open: http://localhost:5000

## Calculators

1. **Z-Factor** — Hall-Yarborough iterative method + gas density & viscosity (Lee-Gonzalez-Eakin)
2. **Gas Compressibility (Cg)** — Isothermal compressibility via dZ/dP numerical differentiation
3. **FVF & OGIP** — Gas formation volume factor Bg + volumetric OGIP
4. **Material Balance** — p/Z plot with regression for OGIP estimation
5. **Pseudo-Critical Properties** — Standing/Sutton correlations + Wichert-Aziz sour gas correction

## Deploy to Render/Railway/Heroku

Add a `Procfile`:
```
web: gunicorn app:app
```
And install gunicorn: `pip install gunicorn`
