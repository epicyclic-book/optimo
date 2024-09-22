#!/usr/bin/python

#############################
# Wolfram Gearbox - Stage 2 #
#############################

import math

##############
# Parameters #
##############

obj_names = (
  'Center distance',
  'Weight',
  'Efficiency',
  'Slide roll ratio',
  'Specific sliding',
  'Almen product',
  'Flash temperature 1',
  'Flash temperature 2',
  'Lambda ratio'
)

obj_types = (
  'min', 'min', 'max',
  'min', 'min', 'min',
  'min', 'min',
  'max'
)

constr_names = (
  'Top land thickness 1',
  'Top land thickness 2',
  'Contact ratio',
  'Bending stress',
  'Contact stress',
  'Involute interference 1',
  'Involute interference 2',
  'Tip interference'
)

var_names = (
  'z1', 'z3',
  'xn1', 'xn3',
  'mn', 'np'
)

var_types = (
  'int', 'int',
  'float', 'float',
  'float', 'int'
)

plot_data = (
  ('Efficiency', 'Center distance', 'Weight'),
  ('Efficiency', 'Slide roll ratio', 'Lambda ratio'),
  ('Efficiency', 'Specific sliding', 'Almen product'),
  ('Efficiency', 'Flash temperature 1', 'Flash temperature 2')
)

config = {
  'num_vars': 6,
  'num_obj': 9,
  'num_constr': 8,
  'bounds_min': [20, 70, 0.0, 0.0, 0.5, 2],
  'bounds_max': [40, 90, 2.0, 2.0, 2.0, 8],
  'ref_data': [32, 81, 0.536, 1.21, 0.8467, 3],

  'num_iter': 100,
  'num_pop': 200,
  'num_archive': 200,

  'obj_names': obj_names,
  'obj_types': obj_types,
  'constr_names': constr_names,
  'var_names': var_names,
  'var_types': var_types,
  'plot_data': plot_data
}

######################
# Problem Definition #
######################

def Problem(v, obj, con):
  # constants
  alpha_n = math.radians(20)
  beta = math.radians(0)
  b = 16
  N1, N3 = 461.54, 31.17
  P = 0.4

  # variables
  z1, z3 = math.ceil(v[0]), math.ceil(v[1])
  x_n1, x_n3 = v[2], v[3]
  m_n, n_w = v[4], math.ceil(v[5])

  # generalized equations
  m_t = m_n / math.cos(beta)
  alpha_t = math.atan2(math.tan(alpha_n), math.cos(beta))

  inv_alpha_t = math.tan(alpha_t) - alpha_t
  inv_alpha_n = math.tan(alpha_n) - alpha_n

  inv_alpha_wt2 = EP2 = 2 * math.tan(alpha_n) * ((x_n3 - x_n1) / (z3 - z1)) + inv_alpha_t

  alpha_wt2 = math.pow(EP2, 1/3) / (0.693357 + 0.192484 * math.pow(EP2, 2/3))

  d1 = z1 * m_n / math.cos(beta)
  d3 = z3 * m_n / math.cos(beta)

  d_b1 = d1 * math.cos(alpha_t)
  d_b3 = d3 * math.cos(alpha_t)

  d_w1 = d_b1 / math.cos(alpha_wt2)
  d_w3 = d_b3 / math.cos(alpha_wt2)

  a = (d3 - d1) / 2
  a_dash = (d_w3 - d_w1) / 2
  del_a = a_dash - a

  y2 = a_dash / m_n - (z3 - z1) / (2 * math.cos(beta))

  h_a1_dash = (1 - y2 + x_n3) * m_n
  h1_dash = (2.25 + y2 - (x_n3 - x_n1)) * m_n

  d_a1_dash = d1 + 2 * h_a1_dash
  d_f1_dash = d_a1_dash - 2 * h1_dash
  d_a3_dash = 2 * (a_dash + d_f1_dash / 2 + m_n / 4)

  r_b1, r_b3 = d_b1 / 2, d_b3 / 2
  r_w1, r_w3 = d_w1 / 2, d_w3 / 2
  r_a1_dash, r_a3_dash = d_a1_dash / 2, d_a3_dash / 2

  alpha_at1 = math.acos(d1 * math.cos(alpha_t) / d_a1_dash)
  alpha_at3 = math.acos(d3 * math.cos(alpha_t) / d_a3_dash)

  # objective 1 - minimize center distance
  obj.append(m_n * (z3 - z1) / 2)

  # objective 2 - minimize weight
  ro = 7.86e-6

  A1 = A3 = b * 2 * m_n
  d1_bar = d_w1 - 2 * m_n
  d3_bar = d_w3 + 2 * m_n

  W1 = ro * math.pi * d1_bar * A1 * n_w
  W3 = ro * math.pi * d3_bar * A3

  obj.append(W1 + W3)

  # objective 3 - maximize efficiency
  beta_b = math.asin(math.sin(beta) * math.cos(alpha_n))

  p_b = math.pi * d_b1 / z1
  p_bt = p_b / math.cos(beta)

  LF = math.sqrt(math.pow(r_a1_dash, 2) - math.pow(r_b1, 2))
  KE = math.sqrt(math.pow(r_a3_dash, 2) - math.pow(r_b3, 2))

  LP = math.sqrt(math.pow(r_w1, 2) - math.pow(r_b1, 2))
  KP = math.sqrt(math.pow(r_w3, 2) - math.pow(r_b3, 2))
  KL = a_dash * math.sin(alpha_wt2)

  EP = KP - KE
  FP = LF - LP
  KF = KP + FP
  LE = KE - KL
  EF = LF - LE

  eps1_dash = FP / p_bt
  eps3 = EP / p_bt

  eps_a2 = eps1_dash + eps3

  i2 = z3 / z1

  Hv2 = (math.pi / (z1 * math.cos(beta_b))) * ((1 + i2) / i2) * \
        (1 - eps_a2 + math.pow(eps1_dash, 2) + math.pow(eps3, 2))

  R_a1 = R_a2 = 0.8
  eta_oil = 50
  d = 0.0651

  ro1_dash = r_w1 * math.sin(alpha_wt2) / math.cos(beta_b)
  ro3 = r_w3 * math.sin(alpha_wt2) / math.cos(beta_b)

  ro_redc2 = ro1_dash * ro3 / ((ro3 - ro1_dash) * math.cos(beta_b))

  T1 = 60000 * P / (2 * math.pi * N1 * n_w)

  F_tb1 = F_tb2 = 1000 * T1 / (r_b1)

  X_L1 = 1 / math.pow(F_tb1 / b, d)
  X_L2 = 1 / math.pow(F_tb2 / b, d)

  omega1 = 2 * math.pi * N1 / 60
  omega2 = 0
  omega3 = 2 * math.pi * N3 / 60

  v_t1 = r_w1 * omega1 / 1000
  v_t2 = 0
  v_t3 = r_w3 * omega3 / 1000

  gamma_epc2 = 2 * v_t1 * math.sin(alpha_wt2)

  mu_m2 = 0.048 * math.pow((F_tb2 / b) / (gamma_epc2 * ro_redc2), 0.2) * \
                  math.pow(eta_oil, -0.05) * math.pow(R_a2, 0.25) * X_L2

  eta = 1 - mu_m2 * Hv2

  obj.append(-(eta * 100))

  # objective 4: minimize slide roll ratio
  v_r1 = ro1_dash * omega1
  v_r3 = ro3 * omega3

  v_s3 = v_r1 - v_r3
  v_s4 = v_r3 - v_r1
  v_e  = v_r1 + v_r3

  obj.append(math.fabs(v_s3) / v_e)

  # objective 5: minimize specific sliding
  gamma_1max_dash = (math.tan(alpha_at1) - math.tan(alpha_wt2)) * (i2 + 1) / \
                    ((i2 + 1) * math.tan(alpha_wt2) - math.tan(alpha_at1))
  gamma_3max      = (math.tan(alpha_at3) - math.tan(alpha_wt2)) * ((i2 + 1) / i2) / \
                    (((i2 + 1) / i2) * math.tan(alpha_wt2) - math.tan(alpha_at3))

  obj.append(math.fabs(gamma_3max - gamma_1max_dash))

  # objective 6: minimize almen product
  v1 = v2 = v3 = 0.3
  E1 = E2 = E3 = 210000
  Z_E1 = Z_E2 = 191.65

  ro_1a_dash = LE / math.cos(beta_b)
  ro_1r_dash = LF / math.cos(beta_b)
  ro_3a = KE / math.cos(beta_b)
  ro_3r = KF / math.cos(beta_b)

  F_na_ddash = 1000 * T1 / (math.sqrt(math.pow(r_b1, 2) + math.pow(LE, 2)) *
                            math.cos(alpha_t) * b)
  F_nr_ddash = 1000 * T1 / (r_a1_dash * math.cos(alpha_t) * b)

  sigma_Ha2 = Z_E2 * math.sqrt(F_na_ddash * (1 / ro_1a_dash - 1 / ro_3a))
  sigma_Hr2 = Z_E2 * math.sqrt(F_nr_ddash * (1 / ro_1r_dash - 1 / ro_3r))

  v_sa2 = EP * (omega2 - omega3) / 1000
  v_sr2 = FP * (omega2 - omega3) / 1000

  obj.append(math.fabs(sigma_Ha2 * v_sa2 - sigma_Hr2 * v_sr2))

  # objective 7-8: minimize flash temperature
  B_m1 = B_m2 = 12.083
  E1_dash = E2_dash = 230770
  K = 0.8

  eps_alpha2_dash = EF / (math.pi * m_t * math.cos(alpha_t))
  eps_beta = b * math.sin(beta) / (math.pi * m_n)

  X_r2 = 1 / eps_alpha2_dash

  ro_redc_2a = ro_1a_dash * ro_3a / ((ro_3a - ro_1a_dash) * math.cos(beta_b))
  ro_redc_2r = ro_1r_dash * ro_3r / ((ro_3r - ro_1r_dash) * math.cos(beta_b))

  b_Ha2 = math.sqrt(4 * F_na_ddash * ro_redc_2a / (math.pi * E2_dash))
  b_Hr2 = math.sqrt(4 * F_nr_ddash * ro_redc_2r / (math.pi * E2_dash))

  v_r1a_ddash = ro_1a_dash * omega1 / 1000
  v_r1r_ddash = ro_1r_dash * omega1 / 1000
  v_r3a_ddash = ro_3a * omega3 / 1000
  v_r3r_ddash = ro_3r * omega3 / 1000

  theta_fla2 = 31.62 * K * mu_m2 * (X_r2 * F_na_ddash / math.sqrt(b_Ha2)) * \
               math.fabs(v_r1a_ddash - v_r3a_ddash) / (B_m1 * math.sqrt(v_r1a_ddash) + \
                                                       B_m2 * math.sqrt(v_r3a_ddash))
  theta_flr2 = 31.62 * K * mu_m2 * (X_r2 * F_nr_ddash / math.sqrt(b_Hr2)) * \
               math.fabs(v_r1r_ddash - v_r3r_ddash) / (B_m1 * math.sqrt(v_r1r_ddash) + \
                                                       B_m2 * math.sqrt(v_r3r_ddash))

  obj.append(theta_fla2)
  obj.append(theta_flr2)

  # objective 9: maximize lambda ratio
  R_a1_dash = R_a2_dash = 0.1
  alpha1_dash = alpha2_dash = 2.19e-8
  eta_o1 = eta_o2 = 90 * 1e-3
  E1_ddash = E2_ddash = 230770 * 1e6

  R_i2 = r_w1 * r_w3 * math.sin(alpha_wt2) / (1000 * (r_w3 - r_w1))

  F_n1 = F_n2 = 1000 * T1 / (r_w1)

  A_ddash = 1000 * F_n2 / (b * E2_ddash * R_i2)
  B_ddash = eta_o2 * (v_t1 - v_t3) * math.cos(alpha_n) / (E2_ddash * R_i2 * 2)
  C_ddash = alpha2_dash * E2_ddash

  h_min2 = 1.714 * R_i2 * math.pow(A_ddash, -0.128) * \
                          math.pow(B_ddash,  0.694) * \
                          math.pow(C_ddash,  0.568)

  obj.append(-h_min2 * 1e6 / R_a2_dash)

  # constraint 1-2: top land thickness
  inv_alpha_at1 = math.tan(alpha_at1) - alpha_at1
  inv_alpha_at3 = math.tan(alpha_at3) - alpha_at3

  S1 = math.pi * m_n / 2 + 2 * x_n1 * m_n * math.tan(alpha_n)
  S3 = math.pi * m_n / 2 - 2 * x_n3 * m_n * math.tan(alpha_n)

  S_a1 = d_a1_dash * (S1 / d1 + inv_alpha_n - inv_alpha_at1)
  S_a3 = d_a3_dash * (S3 / d3 - inv_alpha_n + inv_alpha_at3)

  con.append(S_a1 - 0.2 * m_n)
  con.append(S_a3 - 0.2 * m_n)

  # constraint 3: contact ratio
  con.append(eps_alpha2_dash + eps_beta - 1.2)

  # constraint 4: bending stress
  y_F = 0.84
  y_eps = 1.0 / eps_alpha2_dash
  y_beta = 1.0 - beta / math.radians(120)
  k_v = k_L = k_FX = 1.0
  k_o = 1.75
  S_F = 1.2

  sigma_F = (y_F * y_eps * y_beta / (m_n * b)) * \
            (k_v * k_o / (k_L * k_FX)) * S_F

  sigma_F1 = F_n1 * sigma_F

  con.append(400 - sigma_F1)

  # constraint 5: contact stress
  Z_H2 = math.sqrt(2 * math.cos(beta_b) / math.tan(alpha_wt2)) / math.cos(alpha_t)
  Z_M = 60.62
  Z_beta = 1.0
  K_HL = K_HX = 1.0
  Z_L = Z_W = 1.0
  Z_R = 0.92
  Z_V = 0.96
  K_Hbeta = 1.0
  K_V = 1.0
  K_O = 1.0
  S_H = 1.15

  F_t2 = F_n1 / 9.81

  if eps_beta <= 1:
    Z_eps2 = math.sqrt(1 - eps_beta + eps_beta / eps_alpha2_dash)
  else:
    Z_eps2 = math.sqrt(1 / eps_alpha2_dash)

  sigma_H = (Z_M * Z_beta / (K_HL * Z_L * Z_R * Z_V * Z_W * K_HX)) * \
            math.sqrt(K_Hbeta * K_V * K_O) * S_H * 9.81

  sigma_H2 = math.sqrt((F_t2 * (i2 - 1)) / (d_w1 * b * i2)) * Z_H2 * Z_eps2 * sigma_H

  con.append(2000 - sigma_H2)

  # constraint 6-7: involute interference
  x1 = 1.0 - math.tan(alpha_at3) / math.tan(alpha_wt2)
  x2 = z1 / z3

  con.append(KE - KL)
  con.append(x2 - x1)

  # constraint 8: tip interference
  cos_eps1 = (math.pow(r_a3_dash, 2) - math.pow(r_a1_dash, 2) - \
              math.pow(a_dash, 2)) / (2 * a_dash * r_a1_dash)

  eps1 = math.acos(cos_eps1)
  eps2 = (z1 / z3) * eps1

  delta = (z1 / z3) * (inv_alpha_at1 - inv_alpha_t) + inv_alpha_t
  theta_t1 = math.asin(r_a1_dash * math.sin(eps1) / r_a3_dash) - eps2

  x1 = inv_alpha_at3
  x2 = delta - theta_t1

  con.append(x2 - x1)

