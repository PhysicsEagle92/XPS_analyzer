import numpy as np

Te3d1_peak1_area = np.array([381766.21325268596])
Te3d1_peak2_area = np.array([562329.7210782899])
Te3d1_area = Te3d1_peak1_area + Te3d1_peak2_area

In3p4_grazing = np.array([136241.56583765792])
In3p4_normal = np.array([170796.29389255476])

Te3d1_peak2_area_grazing = np.array([383145.00354032265])
Te3d1_peak1_area_grazing = np.array([533047.7269739407])
Te3d1_grazing_area = Te3d1_peak1_area_grazing + Te3d1_peak2_area_grazing

Te3d1_peak2_area_normal = np.array([379319.2220795632])
Te3d1_peak2_area_normal = np.array([582374.1223294177])
Te3d1_normal_area = Te3d1_peak2_area_normal + Te3d1_peak2_area_normal

Mn3s2_area = np.array([24997.825934221575])

Mn2p2_set1_area = np.array([384911.25350741186])
Mn2p2_set2_area = np.array([35509.179796332915])
Mn2p_area = Mn2p2_set1_area + Mn2p2_set2_area

Mn2p_grazing_set1_area = np.array([147044.47614646782])
Mn2p_grazing_set2_area = np.array([21632.99580676744])
Mn2p_area_grazing = Mn2p_grazing_set1_area + Mn2p_grazing_set2_area

Mn2p_crosssection = np.array([0.1878,0.188,0.1881,1.447,1.447,1.446])
Mn3s_crosssection = np.array([0.009195,0.009145,0.009146])
Te3d_crosssection = np.array([0.4205,0.4205,0.4208,1.185,1.186,1.184])
In3p_crosssection = np.array([0.183,0.1831,0.1835,1.587,1.587,1.578])

print("Mn/Te ratio = ",(Mn2p_area/Mn2p_crosssection[0])/(Te3d1_area/Te3d_crosssection[0]))
print("Mn/Te ratio = ",(Mn3s2_area/Mn3s_crosssection[0])/(Te3d1_area/Te3d_crosssection[0]))

print("Te/In ratio normal = ", (Te3d1_normal_area/Te3d_crosssection[0])/(In3p4_normal/In3p_crosssection[0]))
print("Te/In ratio grazing = ", (Te3d1_grazing_area/Te3d_crosssection[0])/(In3p4_grazing/In3p_crosssection[0]))

print("Mn/Te ratio grazing = ", (Mn2p_area_grazing/Mn2p_crosssection[0])/(Te3d1_grazing_area/Te3d_crosssection[0]))

print("Mn/In ratio grazing = ", (Mn2p_area_grazing/Mn2p_crosssection[0])/(In3p4_grazing/In3p_crosssection[0]))

print("Mn_grazing/Mn ratio = ", (Mn2p_area_grazing/Mn2p_crosssection[0])/(Mn2p_area/Mn2p_crosssection[0]))

print("Te_grazing/Te ratio = ", (Te3d1_grazing_area/Te3d_crosssection[0])/(Te3d1_area/Te3d_crosssection[0]))

print("Mn/In ration = ", (Mn2p_area/Mn2p_crosssection[0])/(In3p4_normal/In3p_crosssection[0]))


