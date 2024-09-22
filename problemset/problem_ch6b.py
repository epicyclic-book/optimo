#!/usr/bin/python

#####################################
# Electric Vehicle Gearbox - Gear 2 #
#####################################

import math

##############
# Parameters #
##############

obj_names = (
  'Center distance',
  'Weight',
  'Efficiency',
  'Specific sliding',
  'Almen product',
  'Flash temperature 1',
  'Flash temperature 2',
  'Lambda ratio'
)

obj_types = (
  'min', 'min', 'max',
  'min', 'min',
  'min', 'min',
  'max'
)

constr_names = (
  'Interference 1',
  'Interference 2',
  'Top land thickness 1',
  'Top land thickness 2',
  'Contact ratio',
  'Bending stress',
  'Contact stress'
)

var_names = (
  'z1', 'z2',
  'xn1', 'xn2',
  'mn', 'np'
)

var_types = (
  'int', 'int',
  'float', 'float',
  'float', 'int'
)

plot_data = (
  ('Efficiency', 'Center distance', 'Weight'),
  ('Efficiency', 'Center distance', 'Lambda ratio'),
  ('Efficiency', 'Specific sliding', 'Almen product'),
  ('Efficiency', 'Flash temperature 1', 'Flash temperature 2')
)

config = {
  'num_vars': 6,
  'num_obj': 8,
  'num_constr': 7,
  'bounds_min': [15, 25, 0.0, 0.0, 1, 2],
  'bounds_max': [25, 40, 0.5, 0.5, 3, 5],
  'ref_data': [19, 32, 0.3, 0.1, 2, 3],

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
  b = 30
  N1, N2 = 1266.8, 0
  P = 100

  # variables
  z1, z2 = math.ceil(v[0]), math.ceil(v[1])
  x_n1, x_n2 = v[2], v[3]
  m_n, n_w = v[4], math.ceil(v[5])

  # generalized equations
  m_t = m_n / math.cos(beta)
  alpha_t = math.atan2(math.tan(alpha_n), math.cos(beta))

  inv_alpha_t = math.tan(alpha_t) - alpha_t
  inv_alpha_n = math.tan(alpha_n) - alpha_n

  inv_alpha_wt1 = EP1 = 2 * math.tan(alpha_n) * ((x_n1 + x_n2) / (z1 + z2)) + inv_alpha_t

  alpha_wt1 = math.pow(EP1, 1/3) / (0.693357 + 0.192484 * math.pow(EP1, 2/3))

  d1 = z1 * m_n / math.cos(beta)
  d2 = z2 * m_n / math.cos(beta)

  d_b1 = d1 * math.cos(alpha_t)
  d_b2 = d2 * math.cos(alpha_t)

  d_w1 = d_b1 / math.cos(alpha_wt1)
  d_w2 = d_b2 / math.cos(alpha_wt1)

  a = (d1 + d2) / 2
  a_dash = (d_w1 + d_w2) / 2
  del_a = a_dash - a

  y1 = a_dash / m_n - (z1 + z2) / (2 * math.cos(beta))

  h_a1_dash = (1 + y1 - x_n2) * m_n
  h_a2_dash = (1 + y1 - x_n1) * m_n

  d_a1_dash = d1 + 2 * h_a1_dash
  d_a2_dash = d2 + 2 * h_a2_dash

  r_b1, r_b2 = d_b1 / 2, d_b2 / 2
  r_w1, r_w2 = d_w1 / 2, d_w2 / 2
  r_a1_dash, r_a2_dash = d_a1_dash / 2, d_a2_dash / 2

  alpha_at1 = math.acos(d1 * math.cos(alpha_t) / d_a1_dash)
  alpha_at2 = math.acos(d2 * math.cos(alpha_t) / d_a2_dash)

  # objective 1 - minimize center distance
  obj.append(m_n * (z1 + z2) / 2)

  # objective 2 - minimize weight
  ro = 7.86e-6

  A1 = A2 = b * 2 * m_n
  d1_bar = d_w1 - 2 * m_n
  d2_bar = d_w2 - 2 * m_n

  W1 = ro * math.pi * d1_bar * A1 * n_w
  W2 = ro * math.pi * d2_bar * A2

  obj.append(W1 + W2)

  # objective 3 - maximize efficiency
  beta_b = math.asin(math.sin(beta) * math.cos(alpha_n))

  p_b = math.pi * d_b1 / z1
  p_bt = p_b / math.cos(beta)

  AP = math.sqrt(math.pow(r_a2_dash, 2) - math.pow(r_b2, 2)) - r_w2 * math.sin(alpha_wt1)
  BP = math.sqrt(math.pow(r_a1_dash, 2) - math.pow(r_b1, 2)) - r_w1 * math.sin(alpha_wt1)
  LF = math.sqrt(math.pow(r_a1_dash, 2) - math.pow(r_b1, 2))
  MP = r_w1 * math.sin(alpha_wt1)
  NP = r_w2 * math.sin(alpha_wt1)

  LP = math.sqrt(math.pow(r_w1, 2) - math.pow(r_b1, 2))

  AM = MP - AP
  FP = LF - LP
  AB = AP + BP

  eps1 = AP / p_bt
  eps2 = BP / p_bt

  eps_a1 = eps1 + eps2

  i1 = z2 / z1

  Hv1 = (math.pi / (z1 * math.cos(beta_b))) * ((1 + i1) / i1) * \
        (1 - eps_a1 + math.pow(eps1, 2) + math.pow(eps2, 2))

  R_a1 = R_a2 = 0.8
  eta_oil = 50
  d = 0.0651

  ro1 = r_w1 * math.sin(alpha_wt1) / math.cos(beta_b)
  ro2 = r_w2 * math.sin(alpha_wt1) / math.cos(beta_b)

  ro_redc1 = ro1 * ro2 / ((ro1 + ro2) * math.cos(beta_b))

  T1 = 60000 * P / (2 * math.pi * N1 * n_w)

  F_tb1 = F_tb2 = 1000 * T1 / (r_b1)

  X_L1 = 1 / math.pow(F_tb1 / b, d)
  X_L2 = 1 / math.pow(F_tb2 / b, d)

  omega1 = 2 * math.pi * N1 / 60
  omega2 = 2 * math.pi * N2 / 60

  v_t1 = r_w1 * omega1 / 1000
  v_t2 = r_w2 * omega2 / 1000

  gamma_epc1 = 2 * v_t1 * math.sin(alpha_wt1)

  mu_m1 = 0.048 * math.pow((F_tb1 / b) / (gamma_epc1 * ro_redc1), 0.2) * \
                  math.pow(eta_oil, -0.05) * math.pow(R_a1, 0.25) * X_L1

  eta = 1 - mu_m1 * Hv1

  obj.append(-(eta * 100))

  # objective 4: minimize specific sliding
  gamma_1max = (math.tan(alpha_at1) - math.tan(alpha_wt1)) * (i1 + 1) / \
               ((i1 + 1) * math.tan(alpha_wt1) - math.tan(alpha_at1))
  gamma_2max = (math.tan(alpha_at2) - math.tan(alpha_wt1)) * ((i1 + 1) / i1) / \
               (((i1 + 1) / i1) * math.tan(alpha_wt1) - math.tan(alpha_at2))

  obj.append(math.fabs(gamma_1max - gamma_2max))

  # objective 5: minimize almen product
  v1 = v2 = v3 = 0.3
  E1 = E2 = E3 = 210000
  Z_E1 = Z_E2 = 191.65

  ro_1a = (MP - AP) / math.cos(beta_b)
  ro_1r = (MP + BP) / math.cos(beta_b)
  ro_2a = (NP + AP) / math.cos(beta_b)
  ro_2r = (NP - BP) / math.cos(beta_b)

  F_na_dash  = 1000 * T1 / (math.sqrt(math.pow(r_b1, 2) + math.pow(AM, 2)) *
                            math.cos(alpha_t) * b)
  F_nr_dash  = 1000 * T1 / (r_a1_dash * math.cos(alpha_t) * b)

  sigma_Ha1 = Z_E1 * math.sqrt(F_na_dash * (1 / ro_1a + 1 / ro_2a))
  sigma_Hr1 = Z_E1 * math.sqrt(F_nr_dash * (1 / ro_1r + 1 / ro_2r))

  v_sa1 = AP * (omega1 + omega2) / 1000
  v_sr1 = BP * (omega1 + omega2) / 1000

  obj.append(math.fabs(sigma_Ha1 * v_sa1 - sigma_Hr1 * v_sr1))

  # objective 6-7: minimize flash temperature
  B_m1 = B_m2 = 12.083
  E1_dash = E2_dash = 230770
  K = 0.8

  eps_alpha1_dash = AB / (math.pi * m_t * math.cos(alpha_t))
  eps_beta = b * math.sin(beta) / (math.pi * m_n)

  X_r1 = 1 / eps_alpha1_dash

  ro_redc_1a = ro_1a * ro_2a / ((ro_1a + ro_2a) * math.cos(beta_b))
  ro_redc_1r = ro_1r * ro_2r / ((ro_1r + ro_2r) * math.cos(beta_b))

  b_Ha1 = math.sqrt(4 * F_na_dash  * ro_redc_1a / (math.pi * E1_dash))
  b_Hr1 = math.sqrt(4 * F_nr_dash  * ro_redc_1r / (math.pi * E1_dash))

  v_r1a_dash = ro_1a * omega1 / 1000
  v_r1r_dash = ro_1r * omega1 / 1000
  v_r2a_dash = ro_2a * omega2 / 1000
  v_r2r_dash = ro_2r * omega2 / 1000

  theta_fla1 = 31.62 * K * mu_m1 * (X_r1 * F_na_dash / math.sqrt(b_Ha1)) * \
               math.fabs(v_r1a_dash - v_r2a_dash) / (B_m1 * math.sqrt(v_r1a_dash) + \
                                                     B_m2 * math.sqrt(v_r2a_dash))
  theta_flr1 = 31.62 * K * mu_m1 * (X_r1 * F_nr_dash / math.sqrt(b_Hr1)) * \
               math.fabs(v_r1r_dash - v_r2r_dash) / (B_m1 * math.sqrt(v_r1r_dash) + \
                                                     B_m2 * math.sqrt(v_r2r_dash))

  obj.append(theta_fla1)
  obj.append(theta_flr1)

  # objective 8: maximize lambda ratio
  R_a1_dash = R_a2_dash = 0.1
  alpha1_dash = alpha2_dash = 2.19e-8
  eta_o1 = eta_o2 = 90 * 1e-3
  E1_ddash = E2_ddash = 230770 * 1e6

  R_i1 = r_w1 * r_w2 * math.sin(alpha_wt1) / (1000 * (r_w1 + r_w2))

  F_n1 = F_n2 = 1000 * T1 / (r_w1)

  A_dash  = 1000 * F_n1 / (b * E1_ddash * R_i1)
  B_dash  = eta_o1 * (v_t1 + v_t2) * math.cos(alpha_n) / (E1_ddash * R_i1 * 2)
  C_dash  = alpha1_dash * E1_ddash

  h_min1 = 1.714 * R_i1 * math.pow(A_dash,  -0.128) * \
                          math.pow(B_dash,   0.694) * \
                          math.pow(C_dash,   0.568)

  obj.append(-h_min1 * 1e6 / R_a1_dash)

  # constraint 1-2: interference
  g_p = math.sqrt(math.pow(r_b1, 2) + math.pow(a_dash * math.sin(alpha_wt1), 2)) - r_a1_dash
  g_s = math.sqrt(math.pow(r_b2, 2) + math.pow(a_dash * math.sin(alpha_wt1), 2)) - r_a2_dash

  con.append(g_p)
  con.append(g_s)

  # constraint 3-4: top land thickness
  inv_alpha_at1 = math.tan(alpha_at1) - alpha_at1
  inv_alpha_at2 = math.tan(alpha_at2) - alpha_at2

  S1 = math.pi * m_n / 2 + 2 * x_n1 * m_n * math.tan(alpha_n)
  S2 = math.pi * m_n / 2 + 2 * x_n2 * m_n * math.tan(alpha_n)

  S_a1 = d_a1_dash * (S1 / d1 + inv_alpha_n - inv_alpha_at1)
  S_a2 = d_a2_dash * (S2 / d2 + inv_alpha_n - inv_alpha_at2)

  con.append(S_a1 - 0.2 * m_n)
  con.append(S_a2 - 0.2 * m_n)

  # constraint 5: contact ratio
  con.append(eps_alpha1_dash + eps_beta - 1.2)

  # constraint 6: bending stress
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

  # constraint 7: contact stress
  Z_H1 = math.sqrt(2 * math.cos(beta_b) / math.tan(alpha_wt1)) / math.cos(alpha_t)
  Z_M = 60.62
  Z_beta = 1.0
  K_HL = K_HX = 1.0
  Z_L = Z_W = 1.0
  Z_R = 0.92
  Z_V = 0.98
  K_Hbeta = 1.0
  K_V = 1.0
  K_O = 1.0
  S_H = 1.15

  F_t1 = F_n1 / 9.81

  if eps_beta <= 1:
    Z_eps1 = math.sqrt(1 - eps_beta + eps_beta / eps_alpha1_dash)
  else:
    Z_eps1 = math.sqrt(1 / eps_alpha1_dash)

  sigma_H = (Z_M * Z_beta / (K_HL * Z_L * Z_R * Z_V * Z_W * K_HX)) * \
            math.sqrt(K_Hbeta * K_V * K_O) * S_H * 9.81

  sigma_H1 = math.sqrt((F_t1 * (i1 + 1)) / (d_w1 * b * i1)) * Z_H1 * Z_eps1 * sigma_H

  con.append(2000 - sigma_H1)

