#!/usr/bin/python

##################################
# Wind Turbine Generator Gearbox #
##################################

import math

##############
# Parameters #
##############

obj_names = (
  'Center distance',
  'Weight',
  'Efficiency',
  'Slide roll ratio',
  'Specific sliding 1',
  'Specific sliding 2',
  'Almen product 1',
  'Almen product 2',
  'Flash temperature 1',
  'Flash temperature 2',
  'Flash temperature 3',
  'Flash temperature 4',
  'Lambda ratio 1',
  'Lambda ratio 2'
)

obj_types = (
  'min', 'min', 'max',
  'min', 'min', 'min',
  'min', 'min',
  'min', 'min', 'min', 'min',
  'max', 'max'
)

constr_names = (
  'Interference 1',
  'Interference 2',
  'Top land thickness 1',
  'Top land thickness 2',
  'Top land thickness 3',
  'Contact ratio 1',
  'Contact ratio 2',
  'Bending stress',
  'Contact stress 1',
  'Contact stress 2',
  'Involute interference 1',
  'Involute interference 2',
  'Tip interference',
  'Proper assembly'
)

var_names = (
  'z1', 'z2', 'z3',
  'xn1', 'xn2', 'xn3',
  'mn', 'np'
)

var_types = (
  'int', 'int', 'int',
  'float', 'float', 'float',
  'float', 'int'
)

plot_data = (
  ('Efficiency', 'Center distance', 'Weight'),
  ('Efficiency', 'Center distance', 'Slide roll ratio'),
  ('Efficiency', 'Specific sliding 1', 'Specific sliding 2'),
  ('Efficiency', 'Almen product 1', 'Almen product 2'),
  ('Efficiency', 'Flash temperature 1', 'Flash temperature 2'),
  ('Efficiency', 'Flash temperature 3', 'Flash temperature 4'),
  ('Efficiency', 'Lambda ratio 1', 'Lambda ratio 2')
)

config = {
  'num_vars': 8,
  'num_obj': 14,
  'num_constr': 14,
  'bounds_min': [30, 15, 85, 0.0, 0.0, 0.0, 12, 2],
  'bounds_max': [40, 25, 95, 0.5, 0.5, 2.0, 18, 6],
  'ref_data': [36, 20, 91, 0.3156, 0.4, 1.6429, 16, 3],

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
  beta = math.radians(8)
  b = 310
  N1, N2 = 22.92, 83.25
  P = 2000

  # variables
  z1, z2, z3 = math.ceil(v[0]), math.ceil(v[1]), math.ceil(v[2])
  x_n1, x_n2, x_n3 = v[3], v[4], v[5]
  m_n, n_w = v[6], math.ceil(v[7])

  # generalized equations
  m_t = m_n / math.cos(beta)
  alpha_t = math.atan2(math.tan(alpha_n), math.cos(beta))

  inv_alpha_t = math.tan(alpha_t) - alpha_t
  inv_alpha_n = math.tan(alpha_n) - alpha_n

  inv_alpha_wt1 = EP1 = 2 * math.tan(alpha_n) * ((x_n1 + x_n2) / (z1 + z2)) + inv_alpha_t
  inv_alpha_wt2 = EP2 = 2 * math.tan(alpha_n) * ((x_n3 - x_n1) / (z3 - z1)) + inv_alpha_t

  alpha_wt1 = math.pow(EP1, 1/3) / (0.693357 + 0.192484 * math.pow(EP1, 2/3))
  alpha_wt2 = math.pow(EP2, 1/3) / (0.693357 + 0.192484 * math.pow(EP2, 2/3))

  d1 = z1 * m_n / math.cos(beta)
  d2 = z2 * m_n / math.cos(beta)
  d3 = z3 * m_n / math.cos(beta)

  d_b1 = d1 * math.cos(alpha_t)
  d_b2 = d2 * math.cos(alpha_t)
  d_b3 = d3 * math.cos(alpha_t)

  d_w1 = d_b1 / math.cos(alpha_wt1)
  d_w2 = d_b2 / math.cos(alpha_wt1)
  d_w3 = d_b3 / math.cos(alpha_wt2)

  a = (d1 + d2) / 2
  a_dash = (d_w1 + d_w2) / 2
  del_a = a_dash - a

  y1 = a_dash / m_n - (z1 + z2) / (2 * math.cos(beta))
  y2 = a_dash / m_n - (z3 - z1) / (2 * math.cos(beta))

  h_a1_dash = (1 + y1 - x_n2) * m_n
  h_a2_dash = (1 + y1 - x_n1) * m_n
  h1_dash = (2.25 + y2 - (x_n3 - x_n1)) * m_n

  d_a1_dash = d1 + 2 * h_a1_dash
  d_a2_dash = d2 + 2 * h_a2_dash
  d_f1_dash = d_a1_dash - 2 * h1_dash
  d_a3_dash = 2 * (a_dash + d_f1_dash / 2 + m_n / 4)

  r_b1, r_b2, r_b3 = d_b1 / 2, d_b2 / 2, d_b3 / 2
  r_w1, r_w2, r_w3 = d_w1 / 2, d_w2 / 2, d_w3 / 2
  r_a1_dash, r_a2_dash, r_a3_dash = d_a1_dash / 2, d_a2_dash / 2, d_a3_dash / 2

  alpha_at1 = math.acos(d1 * math.cos(alpha_t) / d_a1_dash)
  alpha_at2 = math.acos(d2 * math.cos(alpha_t) / d_a2_dash)
  alpha_at3 = math.acos(d3 * math.cos(alpha_t) / d_a3_dash)

  # objective 1 - minimize center distance
  obj.append(m_n * (z1 + z2) / 2)

  # objective 2 - minimize weight
  ro = 7.86e-6

  A1 = A2 = A3 = b * 2 * m_n
  d1_bar = d_w1 - 2 * m_n
  d2_bar = d_w2 - 2 * m_n
  d3_bar = d_w3 + 2 * m_n

  W1 = ro * math.pi * d1_bar * A1 * n_w
  W2 = ro * math.pi * d2_bar * A2
  W3 = ro * math.pi * d3_bar * A3

  obj.append(W1 + W2 + W3)

  # objective 3 - maximize efficiency
  beta_b = math.asin(math.sin(beta) * math.cos(alpha_n))

  p_b = math.pi * d_b1 / z1
  p_bt = p_b / math.cos(beta)

  AP = math.sqrt(math.pow(r_a2_dash, 2) - math.pow(r_b2, 2)) - r_w2 * math.sin(alpha_wt1)
  BP = math.sqrt(math.pow(r_a1_dash, 2) - math.pow(r_b1, 2)) - r_w1 * math.sin(alpha_wt1)
  LF = math.sqrt(math.pow(r_a1_dash, 2) - math.pow(r_b1, 2))
  KE = math.sqrt(math.pow(r_a3_dash, 2) - math.pow(r_b3, 2))
  MP = r_w1 * math.sin(alpha_wt1)
  NP = r_w2 * math.sin(alpha_wt1)

  LP = math.sqrt(math.pow(r_w1, 2) - math.pow(r_b1, 2))
  KP = math.sqrt(math.pow(r_w3, 2) - math.pow(r_b3, 2))
  KL = a_dash * math.sin(alpha_wt2)

  AM = MP - AP
  EP = KP - KE
  FP = LF - LP
  KF = KP + FP
  LE = KE - KL
  AB = AP + BP
  EF = LF - LE

  eps1 = AP / p_bt
  eps2 = BP / p_bt
  eps1_dash = FP / p_bt
  eps3 = EP / p_bt

  eps_a1 = eps1 + eps2
  eps_a2 = eps1_dash + eps3

  i1 = z2 / z1
  i2 = z3 / z1

  Hv1 = (math.pi / (z1 * math.cos(beta_b))) * ((1 + i1) / i1) * \
        (1 - eps_a1 + math.pow(eps1, 2) + math.pow(eps2, 2))
  Hv2 = (math.pi / (z1 * math.cos(beta_b))) * ((1 + i2) / i2) * \
        (1 - eps_a2 + math.pow(eps1_dash, 2) + math.pow(eps3, 2))

  R_a1 = R_a2 = 0.8
  eta_oil = 50
  d = 0.0651

  ro1 = r_w1 * math.sin(alpha_wt1) / math.cos(beta_b)
  ro1_dash = r_w1 * math.sin(alpha_wt2) / math.cos(beta_b)
  ro2 = r_w2 * math.sin(alpha_wt1) / math.cos(beta_b)
  ro3 = r_w3 * math.sin(alpha_wt2) / math.cos(beta_b)

  ro_redc1 = ro1 * ro2 / ((ro1 + ro2) * math.cos(beta_b))
  ro_redc2 = ro1_dash * ro3 / ((ro3 - ro1_dash) * math.cos(beta_b))

  T1 = 60000 * P / (2 * math.pi * N1 * n_w)
  T2 = 60000 * P / (2 * math.pi * N2)

  F_tb1 = F_tb2 = 1000 * T1 / (r_b1)

  X_L1 = 1 / math.pow(F_tb1 / b, d)
  X_L2 = 1 / math.pow(F_tb2 / b, d)

  omega1 = 2 * math.pi * N1 / 60
  omega2 = 2 * math.pi * N2 / 60
  omega3 = 0

  v_t1 = r_w1 * omega1 / 1000
  v_t2 = r_w2 * omega2 / 1000
  v_t3 = 0

  gamma_epc1 = 2 * v_t1 * math.sin(alpha_wt1)
  gamma_epc2 = 2 * v_t1 * math.sin(alpha_wt2)

  mu_m1 = 0.048 * math.pow((F_tb1 / b) / (gamma_epc1 * ro_redc1), 0.2) * \
                  math.pow(eta_oil, -0.05) * math.pow(R_a1, 0.25) * X_L1
  mu_m2 = 0.048 * math.pow((F_tb2 / b) / (gamma_epc2 * ro_redc2), 0.2) * \
                  math.pow(eta_oil, -0.05) * math.pow(R_a2, 0.25) * X_L2

  eta1 = 1 - mu_m1 * Hv1
  eta2 = 1 - mu_m2 * Hv2

  eta = 1 - (mu_m1 * Hv1 + mu_m2 * Hv2)

  obj.append(-(eta * 100))

  # objective 4: minimize slide roll ratio
  v_r1 = ro1 * omega1
  v_r2 = ro2 * omega2

  v_s1 = v_r1 - v_r2
  v_s2 = v_r2 - v_r1
  v_e  = v_r1 + v_r2

  obj.append(math.fabs(v_s1) / v_e)

  # objective 5-6: minimize specific sliding
  gamma_1max = (math.tan(alpha_at1) - math.tan(alpha_wt1)) * (i1 + 1) / \
               ((i1 + 1) * math.tan(alpha_wt1) - math.tan(alpha_at1))
  gamma_2max = (math.tan(alpha_at2) - math.tan(alpha_wt1)) * ((i1 + 1) / i1) / \
               (((i1 + 1) / i1) * math.tan(alpha_wt1) - math.tan(alpha_at2))

  gamma_1max_dash = (math.tan(alpha_at1) - math.tan(alpha_wt2)) * (i2 + 1) / \
                    ((i2 + 1) * math.tan(alpha_wt2) - math.tan(alpha_at1))
  gamma_3max      = (math.tan(alpha_at3) - math.tan(alpha_wt2)) * ((i2 + 1) / i2) / \
                    (((i2 + 1) / i2) * math.tan(alpha_wt2) - math.tan(alpha_at3))

  obj.append(math.fabs(gamma_1max - gamma_2max))
  obj.append(math.fabs(gamma_3max - gamma_1max_dash))

  # objective 7-8: minimize almen product
  v1 = v2 = v3 = 0.3
  E1 = E2 = E3 = 210000
  Z_E1 = Z_E2 = 191.65

  ro_1a = (MP - AP) / math.cos(beta_b)
  ro_1r = (MP + BP) / math.cos(beta_b)
  ro_2a = (NP + AP) / math.cos(beta_b)
  ro_2r = (NP - BP) / math.cos(beta_b)
  ro_1a_dash = LE / math.cos(beta_b)
  ro_1r_dash = LF / math.cos(beta_b)
  ro_3a = KE / math.cos(beta_b)
  ro_3r = KF / math.cos(beta_b)

  F_na_dash  = 1000 * T1 / (math.sqrt(math.pow(r_b1, 2) + math.pow(AM, 2)) *
                            math.cos(alpha_t) * b)
  F_nr_dash  = 1000 * T1 / (r_a1_dash * math.cos(alpha_t) * b)

  F_na_ddash = 1000 * T1 / (math.sqrt(math.pow(r_b1, 2) + math.pow(LE, 2)) *
                            math.cos(alpha_t) * b)
  F_nr_ddash = 1000 * T1 / (r_a1_dash * math.cos(alpha_t) * b)

  sigma_Ha1 = Z_E1 * math.sqrt(F_na_dash * (1 / ro_1a + 1 / ro_2a))
  sigma_Hr1 = Z_E1 * math.sqrt(F_nr_dash * (1 / ro_1r + 1 / ro_2r))
  sigma_Ha2 = Z_E2 * math.sqrt(F_na_ddash * (1 / ro_1a_dash - 1 / ro_3a))
  sigma_Hr2 = Z_E2 * math.sqrt(F_nr_ddash * (1 / ro_1r_dash - 1 / ro_3r))

  v_sa1 = AP * (omega1 + omega2) / 1000
  v_sr1 = BP * (omega1 + omega2) / 1000
  v_sa2 = EP * (omega2 - omega3) / 1000
  v_sr2 = FP * (omega2 - omega3) / 1000

  obj.append(math.fabs(sigma_Ha1 * v_sa1 - sigma_Hr1 * v_sr1))
  obj.append(math.fabs(sigma_Ha2 * v_sa2 - sigma_Hr2 * v_sr2))

  # objective 9-12: minimize flash temperature
  B_m1 = B_m2 = 12.083
  E1_dash = E2_dash = 230770
  K = 0.8

  eps_alpha1_dash = AB / (math.pi * m_t * math.cos(alpha_t))
  eps_alpha2_dash = EF / (math.pi * m_t * math.cos(alpha_t))
  eps_beta = b * math.sin(beta) / (math.pi * m_n)

  X_r1 = 1 / eps_alpha1_dash
  X_r2 = 1 / eps_alpha2_dash

  ro_redc_1a = ro_1a * ro_2a / ((ro_1a + ro_2a) * math.cos(beta_b))
  ro_redc_1r = ro_1r * ro_2r / ((ro_1r + ro_2r) * math.cos(beta_b))
  ro_redc_2a = ro_1a_dash * ro_3a / ((ro_3a - ro_1a_dash) * math.cos(beta_b))
  ro_redc_2r = ro_1r_dash * ro_3r / ((ro_3r - ro_1r_dash) * math.cos(beta_b))

  b_Ha1 = math.sqrt(4 * F_na_dash  * ro_redc_1a / (math.pi * E1_dash))
  b_Hr1 = math.sqrt(4 * F_nr_dash  * ro_redc_1r / (math.pi * E1_dash))
  b_Ha2 = math.sqrt(4 * F_na_ddash * ro_redc_2a / (math.pi * E2_dash))
  b_Hr2 = math.sqrt(4 * F_nr_ddash * ro_redc_2r / (math.pi * E2_dash))

  v_r1a_dash = ro_1a * omega1 / 1000
  v_r1r_dash = ro_1r * omega1 / 1000
  v_r2a_dash = ro_2a * omega2 / 1000
  v_r2r_dash = ro_2r * omega2 / 1000

  v_r1a_ddash = ro_1a_dash * omega1 / 1000
  v_r1r_ddash = ro_1r_dash * omega1 / 1000
  v_r3a_ddash = v_r3r_ddash = 0

  theta_fla1 = 31.62 * K * mu_m1 * (X_r1 * F_na_dash / math.sqrt(b_Ha1)) * \
               math.fabs(v_r1a_dash - v_r2a_dash) / (B_m1 * math.sqrt(v_r1a_dash) + \
                                                     B_m2 * math.sqrt(v_r2a_dash))
  theta_flr1 = 31.62 * K * mu_m1 * (X_r1 * F_nr_dash / math.sqrt(b_Hr1)) * \
               math.fabs(v_r1r_dash - v_r2r_dash) / (B_m1 * math.sqrt(v_r1r_dash) + \
                                                     B_m2 * math.sqrt(v_r2r_dash))

  theta_fla2 = 31.62 * K * mu_m2 * (X_r2 * F_na_ddash / math.sqrt(b_Ha2)) * \
               math.fabs(v_r1a_ddash - v_r3a_ddash) / (B_m1 * math.sqrt(v_r1a_ddash) + \
                                                       B_m2 * math.sqrt(v_r3a_ddash))
  theta_flr2 = 31.62 * K * mu_m2 * (X_r2 * F_nr_ddash / math.sqrt(b_Hr2)) * \
               math.fabs(v_r1r_ddash - v_r3r_ddash) / (B_m1 * math.sqrt(v_r1r_ddash) + \
                                                       B_m2 * math.sqrt(v_r3r_ddash))

  obj.append(theta_fla1)
  obj.append(theta_flr1)
  obj.append(theta_fla2)
  obj.append(theta_flr2)

  # objective 13-14: maximize lambda ratio
  R_a1_dash = R_a2_dash = 0.1
  alpha1_dash = alpha2_dash = 2.19e-8
  eta_o1 = eta_o2 = 90 * 1e-3
  E1_ddash = E2_ddash = 230770 * 1e6

  R_i1 = r_w1 * r_w2 * math.sin(alpha_wt1) / (1000 * (r_w1 + r_w2))
  R_i2 = r_w1 * r_w3 * math.sin(alpha_wt2) / (1000 * (r_w3 - r_w1))

  F_n1 = F_n2 = 1000 * T1 / (r_w1)

  A_dash  = 1000 * F_n1 / (b * E1_ddash * R_i1)
  A_ddash = 1000 * F_n2 / (b * E2_ddash * R_i2)
  B_dash  = eta_o1 * (v_t1 + v_t2) * math.cos(alpha_n) / (E1_ddash * R_i1 * 2)
  B_ddash = eta_o2 * (v_t1 - v_t3) * math.cos(alpha_n) / (E2_ddash * R_i2 * 2)
  C_dash  = alpha1_dash * E1_ddash
  C_ddash = alpha2_dash * E2_ddash

  h_min1 = 1.714 * R_i1 * math.pow(A_dash,  -0.128) * \
                          math.pow(B_dash,   0.694) * \
                          math.pow(C_dash,   0.568)
  h_min2 = 1.714 * R_i2 * math.pow(A_ddash, -0.128) * \
                          math.pow(B_ddash,  0.694) * \
                          math.pow(C_ddash,  0.568)

  obj.append(-h_min1 * 1e6 / R_a1_dash)
  obj.append(-h_min2 * 1e6 / R_a2_dash)

  # constraint 1-2: interference
  g_p = math.sqrt(math.pow(r_b1, 2) + math.pow(a_dash * math.sin(alpha_wt1), 2)) - r_a1_dash
  g_s = math.sqrt(math.pow(r_b2, 2) + math.pow(a_dash * math.sin(alpha_wt1), 2)) - r_a2_dash

  con.append(g_p)
  con.append(g_s)

  # constraint 3-5: top land thickness
  inv_alpha_at1 = math.tan(alpha_at1) - alpha_at1
  inv_alpha_at2 = math.tan(alpha_at2) - alpha_at2
  inv_alpha_at3 = math.tan(alpha_at3) - alpha_at3

  S1 = math.pi * m_n / 2 + 2 * x_n1 * m_n * math.tan(alpha_n)
  S2 = math.pi * m_n / 2 + 2 * x_n2 * m_n * math.tan(alpha_n)
  S3 = math.pi * m_n / 2 - 2 * x_n3 * m_n * math.tan(alpha_n)

  S_a1 = d_a1_dash * (S1 / d1 + inv_alpha_n - inv_alpha_at1)
  S_a2 = d_a2_dash * (S2 / d2 + inv_alpha_n - inv_alpha_at2)
  S_a3 = d_a3_dash * (S3 / d3 - inv_alpha_n + inv_alpha_at3)

  con.append(S_a1 - 0.2 * m_n)
  con.append(S_a2 - 0.2 * m_n)
  con.append(S_a3 - 0.2 * m_n)

  # constraint 6-7: contact ratio
  con.append(eps_alpha1_dash + eps_beta - 1.2)
  con.append(eps_alpha2_dash + eps_beta - 1.2)

  # constraint 8: bending stress
  y_F = 0.84
  y_eps = 1.0 / eps_alpha1_dash
  y_beta = 1.0 - beta / math.radians(120)
  k_v = k_L = k_FX = 1.0
  k_o = 1.75
  S_F = 1.2

  sigma_F = (y_F * y_eps * y_beta / (m_n * b)) * \
            (k_v * k_o / (k_L * k_FX)) * S_F

  sigma_F1 = F_n1 * sigma_F

  con.append(400 - sigma_F1)

  # constraint 9-10: contact stress
  Z_H1 = math.sqrt(2 * math.cos(beta_b) / math.tan(alpha_wt1)) / math.cos(alpha_t)
  Z_H2 = math.sqrt(2 * math.cos(beta_b) / math.tan(alpha_wt2)) / math.cos(alpha_t)
  Z_M = 60.62
  Z_beta = 1.0
  K_HL = K_HX = 1.0
  Z_L = Z_W = 1.0
  Z_R = 0.92
  Z_V = 0.97
  K_Hbeta = 1.0
  K_V = 1.0
  K_O = 1.0
  S_H = 1.15

  F_t1 = F_n1 / 9.81

  if eps_beta <= 1:
    Z_eps1 = math.sqrt(1 - eps_beta + eps_beta / eps_alpha1_dash)
    Z_eps2 = math.sqrt(1 - eps_beta + eps_beta / eps_alpha2_dash)
  else:
    Z_eps1 = math.sqrt(1 / eps_alpha1_dash)
    Z_eps2 = math.sqrt(1 / eps_alpha2_dash)

  sigma_H = (Z_M * Z_beta / (K_HL * Z_L * Z_R * Z_V * Z_W * K_HX)) * \
            math.sqrt(K_Hbeta * K_V * K_O) * S_H * 9.81

  sigma_H1 = math.sqrt((F_t1 * (i1 + 1)) / (d_w1 * b * i1)) * Z_H1 * Z_eps1 * sigma_H
  sigma_H2 = math.sqrt((F_t1 * (i2 - 1)) / (d_w1 * b * i2)) * Z_H2 * Z_eps2 * sigma_H

  con.append(2000 - sigma_H1)
  con.append(2000 - sigma_H2)

  # constraint 11-12: involute interference
  x1 = 1.0 - math.tan(alpha_at3) / math.tan(alpha_wt2)
  x2 = z1 / z3

  con.append(KE - KL)
  con.append(x2 - x1)

  # constraint 13: tip interference
  cos_eps1 = (math.pow(r_a3_dash, 2) - math.pow(r_a1_dash, 2) - \
              math.pow(a_dash, 2)) / (2 * a_dash * r_a1_dash)

  eps1 = math.acos(cos_eps1)
  eps2 = (z1 / z3) * eps1

  delta = (z1 / z3) * (inv_alpha_at1 - inv_alpha_t) + inv_alpha_t
  theta_t1 = math.asin(r_a1_dash * math.sin(eps1) / r_a3_dash) - eps2

  x1 = inv_alpha_at3
  x2 = delta - theta_t1

  con.append(x2 - x1)

  # constraint 14: proper assembly
  I = (z2 + z3) / n_w
  con.append(math.floor(I) - I)

