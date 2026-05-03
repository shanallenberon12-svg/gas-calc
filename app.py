from flask import Flask, render_template, request, jsonify
import numpy as np
import math

app = Flask(__name__)


def compute_z(Ppr, Tpr):
    """Hall-Yarborough Z-factor computation (correct sign convention)"""
    t = 1.0 / Tpr
    c = 0.06125 * t * math.exp(-1.2 * (1 - t) ** 2)
    exp_val = 2.18 + 2.82 * t
    y = 0.001
    for _ in range(200):
        yp = max(y, 1e-12)
        F = -c * Ppr + (y + y**2 + y**3 - y**4) / (1 - y)**3 \
            - (14.76*t - 9.76*t**2 + 4.58*t**3)*y**2 \
            + (90.7*t - 242.2*t**2 + 42.4*t**3)*yp**exp_val
        dF = (1 + 4*y + 4*y**2 - 4*y**3 + y**4) / (1 - y)**4 \
             - 2*(14.76*t - 9.76*t**2 + 4.58*t**3)*y \
             + exp_val*(90.7*t - 242.2*t**2 + 42.4*t**3)*yp**(exp_val - 1)
        if abs(dF) < 1e-15:
            break
        y_new = y - F / dF
        y_new = max(1e-8, min(y_new, 0.99))
        if abs(y_new - y) < 1e-10:
            y = y_new
            break
        y = y_new
    return c * Ppr / y

@app.route('/')
def index():
    return render_template('index.html')

# 1. Gas Z-Factor (Hall-Yarborough method)
@app.route('/api/z_factor', methods=['POST'])
def z_factor():
    try:
        data = request.json
        P = float(data['pressure'])       # psia
        T = float(data['temperature'])    # °R
        Ppc = float(data['ppc'])          # pseudo-critical pressure, psia
        Tpc = float(data['tpc'])          # pseudo-critical temperature, °R

        Ppr = P / Ppc
        Tpr = T / Tpc

        Z = compute_z(Ppr, Tpr)

        # Gas density (lb/ft³) using EOS: rho = PM/(ZRT)
        # M = 28.97 * SG, R = 10.73 psia·ft³/(lb-mol·°R)
        sg = float(data.get('sg', 0.65))
        M = 28.97 * sg
        rho = (P * M) / (Z * 10.73 * T)

        # Gas viscosity (Lee-Gonzalez-Eakin correlation)
        K = ((9.4 + 0.02*M) * T**1.5) / (209 + 19*M + T)
        X = 3.5 + (986/T) + 0.01*M
        Y = 2.4 - 0.2*X
        mu_g = 1e-4 * K * math.exp(X * (rho/62.4)**Y)

        return jsonify({
            'z_factor': round(Z, 5),
            'ppr': round(Ppr, 4),
            'tpr': round(Tpr, 4),
            'density': round(rho, 4),
            'viscosity': round(mu_g, 6),
            'status': 'success'
        })
    except Exception as e:
        return jsonify({'status': 'error', 'message': str(e)})


# 2. Gas Compressibility (Cg)
@app.route('/api/compressibility', methods=['POST'])
def compressibility():
    try:
        data = request.json
        P = float(data['pressure'])
        T = float(data['temperature'])
        Ppc = float(data['ppc'])
        Tpc = float(data['tpc'])

        Ppr = P / Ppc
        Tpr = T / Tpc

        # Z at P
        Z1 = compute_z(Ppr, Tpr)

        # dZ/dP via finite difference
        delta = 0.01 * Ppr
        Z2 = compute_z(Ppr + delta, Tpr)
        dZdPpr = (Z2 - Z1) / delta
        dZdP = dZdPpr / Ppc

        Cg = (1/P) - (1/Z1) * dZdP

        return jsonify({
            'cg': f"{Cg:.6e}",
            'cg_value': Cg,
            'z_factor': round(Z1, 5),
            'ppr': round(Ppr, 4),
            'tpr': round(Tpr, 4),
            'status': 'success'
        })
    except Exception as e:
        return jsonify({'status': 'error', 'message': str(e)})


# 3. Gas FVF & Volumetric Reserve
@app.route('/api/fvf', methods=['POST'])
def fvf():
    try:
        data = request.json
        P = float(data['pressure'])
        T = float(data['temperature'])
        Ppc = float(data['ppc'])
        Tpc = float(data['tpc'])
        area = float(data.get('area', 0))         # acres
        thickness = float(data.get('thickness', 0))  # ft
        porosity = float(data.get('porosity', 0))    # fraction
        sw = float(data.get('sw', 0))                # water saturation

        Ppr = P / Ppc
        Tpr = T / Tpc
        Z = compute_z(Ppr, Tpr)

        Bg = 0.02829 * Z * T / P  # res ft³/scf
        Bg_bbl = Bg / 5.615        # res bbl/scf

        # OGIP if reservoir params provided
        ogip = None
        if area > 0 and thickness > 0 and porosity > 0:
            Vb = area * 43560 * thickness  # bulk volume ft³
            Vp = Vb * porosity
            OGIP = Vp * (1 - sw) / Bg      # scf
            ogip = round(OGIP / 1e9, 4)    # Bscf

        return jsonify({
            'bg_cuft': round(Bg, 6),
            'bg_bbl': round(Bg_bbl, 8),
            'z_factor': round(Z, 5),
            'ogip_bscf': ogip,
            'status': 'success'
        })
    except Exception as e:
        return jsonify({'status': 'error', 'message': str(e)})


# 4. Gas Material Balance (p/Z plot)
@app.route('/api/material_balance', methods=['POST'])
def material_balance():
    try:
        data = request.json
        Pi = float(data['pi'])
        Ti = float(data['temperature'])
        Ppc = float(data['ppc'])
        Tpc = float(data['tpc'])
        Gp_list = list(map(float, data['gp'].split(',')))
        P_list = list(map(float, data['pressures'].split(',')))

        Tpr = Ti / Tpc
        pz_points = []

        for P, Gp in zip(P_list, Gp_list):
            Ppr = P / Ppc
            Z = compute_z(Ppr, Tpr)
            pz_points.append({'gp': Gp, 'pz': round(P / Z, 4), 'z': round(Z, 5), 'p': P})

        # Initial p/Z
        Ppr_i = Pi / Ppc
        Zi = compute_z(Ppr_i, Tpr)
        pzi = Pi / Zi

        # Linear regression for OGIP
        x = [pt['gp'] for pt in pz_points]
        y = [pt['pz'] for pt in pz_points]

        if len(x) >= 2:
            n = len(x)
            sx, sy = sum(x), sum(y)
            sxy = sum(xi*yi for xi, yi in zip(x, y))
            sxx = sum(xi**2 for xi in x)
            slope = (n*sxy - sx*sy) / (n*sxx - sx**2)
            intercept = (sy - slope*sx) / n
            OGIP = -intercept / slope if slope != 0 else None
        else:
            OGIP = None
            slope = None
            intercept = pzi

        return jsonify({
            'pz_points': pz_points,
            'pzi': round(pzi, 4),
            'zi': round(Zi, 5),
            'ogip_mmscf': round(OGIP/1e6, 2) if OGIP else None,
            'slope': slope,
            'intercept': round(intercept, 4) if intercept else None,
            'status': 'success'
        })
    except Exception as e:
        return jsonify({'status': 'error', 'message': str(e)})


# 5. Pseudo-Critical Properties from Gas Gravity
@app.route('/api/pseudocritical', methods=['POST'])
def pseudocritical():
    try:
        data = request.json
        sg = float(data['sg'])
        co2 = float(data.get('co2', 0)) / 100
        h2s = float(data.get('h2s', 0)) / 100
        n2 = float(data.get('n2', 0)) / 100

        # Kay's mixing rules (Standing correlation for gas condensate)
        method = data.get('method', 'standing')

        if method == 'sutton':
            Ppc = 756.8 - 131.0*sg - 3.6*sg**2
            Tpc = 169.2 + 349.5*sg - 74.0*sg**2
        else:  # Standing
            Ppc = 677 + 15.0*sg - 37.5*sg**2
            Tpc = 168 + 325*sg - 12.5*sg**2

        # Wichert-Aziz correction for sour gas
        A = h2s + co2
        B = h2s
        if A > 0:
            eps = 120*(A**0.9 - A**1.6) + 15*(B**0.5 - B**4.0)
            Tpc_corr = Tpc - eps
            Ppc_corr = (Ppc * Tpc_corr) / (Tpc + B*(1 - B)*eps)
        else:
            Tpc_corr = Tpc
            Ppc_corr = Ppc

        MW = 28.97 * sg

        return jsonify({
            'ppc': round(Ppc_corr, 2),
            'tpc': round(Tpc_corr, 2),
            'ppc_raw': round(Ppc, 2),
            'tpc_raw': round(Tpc, 2),
            'mw': round(MW, 3),
            'method': method,
            'wichert_applied': A > 0,
            'status': 'success'
        })
    except Exception as e:
        return jsonify({'status': 'error', 'message': str(e)})


if __name__ == '__main__':
    app.run(debug=True, port=5000)
