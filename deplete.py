# @title HTR-PM

import openmc
import openmc.model
import openmc.deplete

import numpy as np
import time, datetime, os, sys, subprocess

from build import (
    region,
    lower_the_reactor,
    seed,
    coretemp,
    fuel,
    enrichment,
    test,
    fuel_volume,
    max_fuel,
    coolant_density,
    mode,
)

wholemattemp = 273 + 20  # Kelvin
np.random.seed(seed)
startstep = 0
regression = 4
cwd = f"{os.getcwd()}/depletion"
try:
    if sys.argv[1] and sys.argv[2]:
        startstep = int(sys.argv[1])
        timestep = int(float(sys.argv[2]))
    else:
        print(sys.argv[2])
        raise ValueError("Tentukan nilai startstep, timestep yang benar")
except Exception:
    raise

gen_time = time.perf_counter()
model_new = None
if startstep != 0:
    model = openmc.Model.from_xml(
        f"{cwd}/step_{startstep-1}/geometry.xml",
        f"{cwd}/step_{startstep-1}/materials.xml",
        f"{cwd}/step_{startstep-1}/settings.xml",
        f"{cwd}/step_{startstep-1}/tallies.xml",
    ) # Take the depleted material and certain cell from this

    # model_new = openmc.Model.from_xml() # Use the new geometry
    if model_new is not None:
        print("=======================================================")
        print("=======================================================")
        print("=======================================================")
        print("================REASIGNING IDS=========================")
        print("=======================================================")
        print("=======================================================")
        print("=======================================================")
        geom_new = model_new.geometry.clone()
        model_new = openmc.Model(
            geom_new, 
            openmc.Materials(list(geom_new.get_all_materials().values())),
            model.settings, model.tallies
        )

else:
    model = openmc.Model.from_xml()
###============================================== ###
###=================== Material ================= ###
###============================================== ###
if fuel == "UCO":
    if enrichment == 'low':
        uco_enrichment = 0.086
        m_fuel_high = openmc.Material(name="Fuel-UCO-Low")
        m_fuel_high.set_density("g/cm3", 10.4)
        m_fuel_high.add_nuclide("U235", 0.885 * uco_enrichment, "wo")
        m_fuel_high.add_nuclide("U238", 0.885 * (1 - uco_enrichment), "wo")
    else:
        uco_enrichment = 0.155
        m_fuel_high = openmc.Material(name="Fuel-UCO-High")
        m_fuel_high.set_density("g/cm3", 10.4)
        m_fuel_high.add_nuclide("U235", 0.885 * uco_enrichment, "wo")
        m_fuel_high.add_nuclide("U238", 0.885 * (1 - uco_enrichment), "wo")
    
    m_fuel_high.add_element("C", 0.019136613, "wo")
    m_fuel_high.add_element("O", 0.095859387, "wo")
    m_fuel_high.add_element("B", 0.000004, "wo")

    m_fuel_high.temperature = coretemp
    m_fuel_high.depletable = True
    m_fuel_high.volume = fuel_volume
    
else:
    if enrichment == "high":
        m_fuel_high = openmc.Material(name="Fuel-UO2-High")
        # Enrichment 15,5%
        m_fuel_high.add_nuclide("U235", 0.136605939, "wo")
        m_fuel_high.add_nuclide("U238", 0.744722699, "wo")
        m_fuel_high.add_element("O",    0.118667362, "wo")

    else:
        m_fuel_high = openmc.Material(name="Fuel-UO2-Low")
        # Enrichment 8,6%
        m_fuel_high.add_nuclide('U235', 0.075802111, 'wo')
        m_fuel_high.add_nuclide('U238', 0.805617787, 'wo')
        m_fuel_high.add_element('O',    0.118576102, 'wo')
        # Reference
        # m_fuel_high.add_nuclide("U234", 1.09e-7, "ao")
        # m_fuel_high.add_nuclide("U235", 1.36e-5, "ao")
        # m_fuel_high.add_nuclide("U238", 1.42e-4, "ao")
        # m_fuel_high.add_nuclide("O16", 3.11e-4, "ao")
        # m_fuel_high.add_nuclide("O17", 1.18e-7, "ao")

    m_fuel_high.add_element("B", 0.000004, "wo")
    m_fuel_high.add_s_alpha_beta("c_U_in_UO2")
    m_fuel_high.set_density("g/cm3", 10.4)
    m_fuel_high.temperature = coretemp
    m_fuel_high.depletable = True
    m_fuel_high.volume = fuel_volume

# Buffer
m_Buffer = openmc.Material(name="Buffer")
m_Buffer.set_density("g/cm3", 1.1)
m_Buffer.add_element('C', 1-0.795E-6, 'wo')
m_Buffer.add_element('B', 0.795E-6, 'wo')
m_Buffer.add_s_alpha_beta("c_Graphite")
m_Buffer.temperature = coretemp

m_PyC = openmc.Material(name="m_PyC")
m_PyC.set_density("g/cm3", 1.9)
m_PyC.add_element('C', 1-0.795E-6, 'wo')
m_PyC.add_element('B', 0.795E-6, 'wo')
m_PyC.add_s_alpha_beta("c_Graphite")
m_PyC.temperature = coretemp

m_SiC = openmc.Material(name="m_SiC")
m_SiC.add_element("C", 12/(12+27.9769265344)-0.795E-6, "wo")
m_SiC.add_element("Si", 27.9769265344/(12+27.9769265344)-0.795E-6, "wo")
m_SiC.add_element("B", 0.795E-6, "wo")
m_SiC.add_s_alpha_beta("c_Graphite")
m_SiC.set_density("g/cm3", 3.18)
m_SiC.temperature = coretemp

m_graphite = openmc.Material(name="m_Graphite_Moderator")
m_graphite.set_density("g/cm3", 1.73)
m_graphite.add_element('C', 1-0.795E-6, 'wo')
m_graphite.add_element('B', 0.795E-6, 'wo')
m_graphite.add_s_alpha_beta("c_Graphite")
m_graphite.temperature = coretemp

m_graphite_pebble = openmc.Material(name="m_Graphite_Pebble")
m_graphite_pebble.set_density("g/cm3", 1.73)
m_graphite_pebble.add_element('C', 1-1e-6, 'wo')
m_graphite_pebble.add_element('B', 1e-6, 'wo')
m_graphite_pebble.add_s_alpha_beta("c_Graphite")
m_graphite_pebble.temperature = coretemp

m_pendingin = openmc.Material(name="m_Pendingin")
if mode == "deplete":
    m_pendingin.add_nuclide("He3", 0.000002)
    m_pendingin.add_nuclide("He4", 0.999998)
    m_pendingin.set_density("g/cm3", coolant_density(coretemp))
else:
    m_pendingin.add_element('H', 2.02935e-6, 'ao')
    m_pendingin.add_element('N', 3.66943e-5, 'ao')
    m_pendingin.add_element('O', 1.08939e-5, 'ao')
    m_pendingin.set_density('g/cm3', 1.146e-3)
m_pendingin.temperature = coretemp


###============================================== ###
###=================== Surface ================== ###
###============================================== ###
h_region = (605 + 562.33) / region
s_reg = [
    openmc.ZPlane(326.17 - lower_the_reactor + (i * h_region))
    for i in range(0, region + 1)
]
if enrichment == "high":
    triso_pitch = 0.109392022682618
    # UCO TRISO
    kernel_r = openmc.Sphere(r=212.5e-4)
    buffer_r = openmc.Sphere(r=(212.5 + 100) * 1e-4)
    IPyC_r = openmc.Sphere(r=(212.5 + 100 + 40) * 1e-4)
    SiC_r = openmc.Sphere(r=(212.5 + 100 + 40 + 35) * 1e-4)
    OPyC_r = openmc.Sphere(r=(212.5 + 100 + 40 + 35 + 40) * 1e-4) # 427.5
else:
    # UO2 TRISO
    triso_pitch = 0.177659234
    kernel_r = openmc.Sphere(r=250e-4)
    buffer_r = openmc.Sphere(r=(250 + 95) * 1e-4)
    IPyC_r = openmc.Sphere(r=(250 + 95 + 40) * 1e-4)
    SiC_r = openmc.Sphere(r=(250 + 95 + 40 + 35) * 1e-4)
    OPyC_r = openmc.Sphere(r=(250 + 95 + 40 + 35 + 40) * 1e-4) # 460
    r_triso = (250 + 95 + 40 + 35 + 40) * 1e-4
pebbles_pitch = 9.05164493902406

graphite_matrix = openmc.Universe(
    name="Graphite-Matrix", cells=[openmc.Cell(name="Graphite-Matrix", fill=m_graphite)]
)
all_coolant = openmc.Universe(
    name="All-Coolant-Universe",
    cells=[openmc.Cell(name="All-Coolant-Cell", fill=m_pendingin)],
)

fcc_pebbles_pitch = pebbles_pitch / 2
top_fcc_pebbles = openmc.ZPlane(fcc_pebbles_pitch)
bottom_fcc_pebbles = openmc.ZPlane(-fcc_pebbles_pitch)
right_fcc_pebbles = openmc.XPlane(fcc_pebbles_pitch)
left_fcc_pebbles = openmc.XPlane(-fcc_pebbles_pitch)
front_fcc_pebbles = openmc.YPlane(fcc_pebbles_pitch)
behind_fcc_pebbles = openmc.YPlane(-fcc_pebbles_pitch)
box_fcc_pebbles = (
    -top_fcc_pebbles
    & +bottom_fcc_pebbles
    & -right_fcc_pebbles
    & +left_fcc_pebbles
    & -front_fcc_pebbles
    & +behind_fcc_pebbles
)

fcc_p_1_coordinates = [-fcc_pebbles_pitch, fcc_pebbles_pitch, fcc_pebbles_pitch]
fcc_p_2_coordinates = [fcc_pebbles_pitch, fcc_pebbles_pitch, fcc_pebbles_pitch]
fcc_p_3_coordinates = [0, 0, fcc_pebbles_pitch]
fcc_p_4_coordinates = [-fcc_pebbles_pitch, -fcc_pebbles_pitch, fcc_pebbles_pitch]
fcc_p_5_coordinates = [fcc_pebbles_pitch, -fcc_pebbles_pitch, fcc_pebbles_pitch]

fcc_p_6_coordinates = [0, fcc_pebbles_pitch, 0]
fcc_p_7_coordinates = [-fcc_pebbles_pitch, 0, 0]
fcc_p_8_coordinates = [fcc_pebbles_pitch, 0, 0]
fcc_p_9_coordinates = [0, -fcc_pebbles_pitch, 0]

fcc_p_10_coordinates = [-fcc_pebbles_pitch, fcc_pebbles_pitch, -fcc_pebbles_pitch]
fcc_p_11_coordinates = [fcc_pebbles_pitch, fcc_pebbles_pitch, -fcc_pebbles_pitch]
fcc_p_12_coordinates = [0, 0, -fcc_pebbles_pitch]
fcc_p_13_coordinates = [-fcc_pebbles_pitch, -fcc_pebbles_pitch, -fcc_pebbles_pitch]
fcc_p_14_coordinates = [fcc_pebbles_pitch, -fcc_pebbles_pitch, -fcc_pebbles_pitch]

s_fcc_f_1 = openmc.Sphere(
    fcc_p_1_coordinates[0], fcc_p_1_coordinates[1], fcc_p_1_coordinates[2], r=2.5
)
s_fcc_f_2 = openmc.Sphere(
    fcc_p_2_coordinates[0], fcc_p_2_coordinates[1], fcc_p_2_coordinates[2], r=2.5
)
s_fcc_f_3 = openmc.Sphere(
    fcc_p_3_coordinates[0], fcc_p_3_coordinates[1], fcc_p_3_coordinates[2], r=2.5
)
s_fcc_f_4 = openmc.Sphere(
    fcc_p_4_coordinates[0], fcc_p_4_coordinates[1], fcc_p_4_coordinates[2], r=2.5
)
s_fcc_f_5 = openmc.Sphere(
    fcc_p_5_coordinates[0], fcc_p_5_coordinates[1], fcc_p_5_coordinates[2], r=2.5
)
s_fcc_f_6 = openmc.Sphere(
    fcc_p_6_coordinates[0], fcc_p_6_coordinates[1], fcc_p_6_coordinates[2], r=2.5
)
s_fcc_f_7 = openmc.Sphere(
    fcc_p_7_coordinates[0], fcc_p_7_coordinates[1], fcc_p_7_coordinates[2], r=2.5
)
s_fcc_f_8 = openmc.Sphere(
    fcc_p_8_coordinates[0], fcc_p_8_coordinates[1], fcc_p_8_coordinates[2], r=2.5
)
s_fcc_f_9 = openmc.Sphere(
    fcc_p_9_coordinates[0], fcc_p_9_coordinates[1], fcc_p_9_coordinates[2], r=2.5
)
s_fcc_f_10 = openmc.Sphere(
    fcc_p_10_coordinates[0], fcc_p_10_coordinates[1], fcc_p_10_coordinates[2], r=2.5
)
s_fcc_f_11 = openmc.Sphere(
    fcc_p_11_coordinates[0], fcc_p_11_coordinates[1], fcc_p_11_coordinates[2], r=2.5
)
s_fcc_f_12 = openmc.Sphere(
    fcc_p_12_coordinates[0], fcc_p_12_coordinates[1], fcc_p_12_coordinates[2], r=2.5
)
s_fcc_f_13 = openmc.Sphere(
    fcc_p_13_coordinates[0], fcc_p_13_coordinates[1], fcc_p_13_coordinates[2], r=2.5
)
s_fcc_f_14 = openmc.Sphere(
    fcc_p_14_coordinates[0], fcc_p_14_coordinates[1], fcc_p_14_coordinates[2], r=2.5
)

s_fcc_g_1 = openmc.Sphere(
    fcc_p_1_coordinates[0], fcc_p_1_coordinates[1], fcc_p_1_coordinates[2], r=3
)
s_fcc_g_2 = openmc.Sphere(
    fcc_p_2_coordinates[0], fcc_p_2_coordinates[1], fcc_p_2_coordinates[2], r=3
)
s_fcc_g_3 = openmc.Sphere(
    fcc_p_3_coordinates[0], fcc_p_3_coordinates[1], fcc_p_3_coordinates[2], r=3
)
s_fcc_g_4 = openmc.Sphere(
    fcc_p_4_coordinates[0], fcc_p_4_coordinates[1], fcc_p_4_coordinates[2], r=3
)
s_fcc_g_5 = openmc.Sphere(
    fcc_p_5_coordinates[0], fcc_p_5_coordinates[1], fcc_p_5_coordinates[2], r=3
)
s_fcc_g_6 = openmc.Sphere(
    fcc_p_6_coordinates[0], fcc_p_6_coordinates[1], fcc_p_6_coordinates[2], r=3
)
s_fcc_g_7 = openmc.Sphere(
    fcc_p_7_coordinates[0], fcc_p_7_coordinates[1], fcc_p_7_coordinates[2], r=3
)
s_fcc_g_8 = openmc.Sphere(
    fcc_p_8_coordinates[0], fcc_p_8_coordinates[1], fcc_p_8_coordinates[2], r=3
)
s_fcc_g_9 = openmc.Sphere(
    fcc_p_9_coordinates[0], fcc_p_9_coordinates[1], fcc_p_9_coordinates[2], r=3
)
s_fcc_g_10 = openmc.Sphere(
    fcc_p_10_coordinates[0], fcc_p_10_coordinates[1], fcc_p_10_coordinates[2], r=3
)
s_fcc_g_11 = openmc.Sphere(
    fcc_p_11_coordinates[0], fcc_p_11_coordinates[1], fcc_p_11_coordinates[2], r=3
)
s_fcc_g_12 = openmc.Sphere(
    fcc_p_12_coordinates[0], fcc_p_12_coordinates[1], fcc_p_12_coordinates[2], r=3
)
s_fcc_g_13 = openmc.Sphere(
    fcc_p_13_coordinates[0], fcc_p_13_coordinates[1], fcc_p_13_coordinates[2], r=3
)
s_fcc_g_14 = openmc.Sphere(
    fcc_p_14_coordinates[0], fcc_p_14_coordinates[1], fcc_p_14_coordinates[2], r=3
)

# Graphite Pebbles
fcc_g_1 = openmc.Cell(name="FCC-G-1", fill=m_graphite_pebble, region=-s_fcc_g_1)
fcc_g_2 = openmc.Cell(name="FCC-G-2", fill=m_graphite_pebble, region=-s_fcc_g_2)
fcc_g_3 = openmc.Cell(name="FCC-G-3", fill=m_graphite_pebble, region=-s_fcc_g_3)
fcc_g_4 = openmc.Cell(name="FCC-G-4", fill=m_graphite_pebble, region=-s_fcc_g_4)
fcc_g_5 = openmc.Cell(name="FCC-G-5", fill=m_graphite_pebble, region=-s_fcc_g_5)
fcc_g_6 = openmc.Cell(name="FCC-G-6", fill=m_graphite_pebble, region=-s_fcc_g_6)
fcc_g_7 = openmc.Cell(name="FCC-G-7", fill=m_graphite_pebble, region=-s_fcc_g_7)
fcc_g_8 = openmc.Cell(name="FCC-G-8", fill=m_graphite_pebble, region=-s_fcc_g_8)
fcc_g_9 = openmc.Cell(name="FCC-G-9", fill=m_graphite_pebble, region=-s_fcc_g_9)
fcc_g_10 = openmc.Cell(name="FCC-G-10", fill=m_graphite_pebble, region=-s_fcc_g_10)
fcc_g_11 = openmc.Cell(name="FCC-G-11", fill=m_graphite_pebble, region=-s_fcc_g_11)
fcc_g_12 = openmc.Cell(name="FCC-G-12", fill=m_graphite_pebble, region=-s_fcc_g_12)
fcc_g_13 = openmc.Cell(name="FCC-G-13", fill=m_graphite_pebble, region=-s_fcc_g_13)
fcc_g_14 = openmc.Cell(name="FCC-G-14", fill=m_graphite_pebble, region=-s_fcc_g_14)

fcc_he_inter1 = openmc.Cell(
    name="Helium-Inter",
    fill=all_coolant,
    region=+s_fcc_g_1
    & +s_fcc_g_2
    & +s_fcc_g_3
    & +s_fcc_g_4
    & +s_fcc_g_5
    & +s_fcc_g_6
    & +s_fcc_g_7
    & +s_fcc_g_8
    & +s_fcc_g_9
    & +s_fcc_g_10
    & +s_fcc_g_11
    & +s_fcc_g_12
    & +s_fcc_g_13
    & +s_fcc_g_14,
)

fcc_g_u_ = openmc.Universe(
    name="FCC-G-Raw",
    cells=[
        fcc_g_1,
        fcc_g_2,
        fcc_g_3,
        fcc_g_4,
        fcc_g_5,
        fcc_g_6,
        fcc_g_7,
        fcc_g_8,
        fcc_g_9,
        fcc_g_10,
        fcc_g_11,
        fcc_g_12,
        fcc_g_13,
        fcc_g_14,
        fcc_he_inter1,
    ],
)
fcc_g_c = openmc.Cell(
    name="FCC-Graphite-Pebbles", fill=fcc_g_u_, region=box_fcc_pebbles
)
fcc_g_u = openmc.Universe(name="FCC-Graphite-Pebbles", cells=[fcc_g_c])

graphite_pebbles_lattice = openmc.RectLattice(name="FCC-Graphite-Pebbles")
graphite_pebbles_lattice.lower_left = (-160, -160, 325 - lower_the_reactor)
graphite_pebbles_lattice.pitch = (pebbles_pitch, pebbles_pitch, pebbles_pitch)
graphite_pebbles_lattice.outer = all_coolant
graphite_pebbles_lattice.universes = [[[fcc_g_u] * 50] * 50] * 70


###============================================== ###
###================= Fuel Pebble ================ ###
###============================================== ###
def build_fuel_lattice(f_mat: openmc.Material, i: int = 1):
    triso_particle = openmc.Universe(
        name="Triso-Particle-High",
        cells=[
            openmc.Cell(name=f"Kernel-High-{i}", fill=f_mat, region=-kernel_r),
            openmc.Cell(name="Buffer", fill=m_Buffer, region=+kernel_r & -buffer_r),
            openmc.Cell(name="IPyC", fill=m_PyC, region=+buffer_r & -IPyC_r),
            openmc.Cell(name="m_SiC", fill=m_SiC, region=+IPyC_r & -SiC_r),
            openmc.Cell(name="OPyC", fill=m_PyC, region=+SiC_r & -OPyC_r),
            openmc.Cell(name="All-Graphite", fill=m_graphite, region=+OPyC_r),
        ],
    )

    n_triso = 6
    triso_6x6_lattice = openmc.RectLattice(name="One-Triso-6x6")
    triso_6x6_lattice.lower_left = (
        -n_triso / 2 * triso_pitch,
        -n_triso / 2 * triso_pitch,
        -n_triso / 2 * triso_pitch,
    )
    triso_6x6_lattice.pitch = (triso_pitch, triso_pitch, triso_pitch)
    triso_6x6_lattice.outer = graphite_matrix
    triso_6x6_lattice.universes = [[[triso_particle] * n_triso] * n_triso] * n_triso

    triso_6x6_cell = openmc.Cell(name="TRISO-6x6-Cell", fill=triso_6x6_lattice)
    triso_6x6_universe = openmc.Universe(
        name="TRISO-6x6-Universe", cells=[triso_6x6_cell]
    )

    triso_lattice = openmc.RectLattice(name="One-Triso")
    triso_lattice.lower_left = (
        -fcc_pebbles_pitch,
        -fcc_pebbles_pitch,
        -fcc_pebbles_pitch,
    )
    triso_lattice.pitch = (
        triso_pitch * n_triso,
        triso_pitch * n_triso,
        triso_pitch * n_triso,
    )
    triso_lattice.outer = graphite_matrix
    triso_lattice.universes = [
        [[triso_6x6_universe] * (int(pebbles_pitch / (triso_pitch * n_triso)) + 1)]
        * (int(pebbles_pitch / (triso_pitch * n_triso)) + 1)
    ] * (int(pebbles_pitch / (triso_pitch * n_triso)) + 1)

    # FCC FUEL PEBBLE
    fcc_fi_1_high = openmc.Cell(name="FCC-FI-1", fill=triso_lattice, region=-s_fcc_f_1)
    fcc_fi_2_high = openmc.Cell(name="FCC-FI-2", fill=triso_lattice, region=-s_fcc_f_2)
    fcc_fi_3_high = openmc.Cell(name="FCC-FI-3", fill=triso_lattice, region=-s_fcc_f_3)
    fcc_fi_4_high = openmc.Cell(name="FCC-FI-4", fill=triso_lattice, region=-s_fcc_f_4)
    fcc_fi_5_high = openmc.Cell(name="FCC-FI-5", fill=triso_lattice, region=-s_fcc_f_5)
    fcc_fi_6_high = openmc.Cell(name="FCC-FI-6", fill=triso_lattice, region=-s_fcc_f_6)
    fcc_fi_7_high = openmc.Cell(name="FCC-FI-7", fill=triso_lattice, region=-s_fcc_f_7)
    fcc_fi_8_high = openmc.Cell(name="FCC-FI-8", fill=triso_lattice, region=-s_fcc_f_8)
    fcc_fi_9_high = openmc.Cell(name="FCC-FI-9", fill=triso_lattice, region=-s_fcc_f_9)
    fcc_fi_10_high = openmc.Cell(
        name="FCC-FI-10", fill=triso_lattice, region=-s_fcc_f_10
    )
    fcc_fi_11_high = openmc.Cell(
        name="FCC-FI-11", fill=triso_lattice, region=-s_fcc_f_11
    )
    fcc_fi_12_high = openmc.Cell(
        name="FCC-FI-12", fill=triso_lattice, region=-s_fcc_f_12
    )
    fcc_fi_13_high = openmc.Cell(
        name="FCC-FI-13", fill=triso_lattice, region=-s_fcc_f_13
    )
    fcc_fi_14_high = openmc.Cell(
        name="FCC-FI-14", fill=triso_lattice, region=-s_fcc_f_14
    )

    # Outer Part of Pebble (Graphite)
    fcc_fo_1 = openmc.Cell(
        name="FCC-FO-1", fill=m_graphite, region=-s_fcc_g_1 & +s_fcc_f_1
    )
    fcc_fo_2 = openmc.Cell(
        name="FCC-FO-2", fill=m_graphite, region=-s_fcc_g_2 & +s_fcc_f_2
    )
    fcc_fo_3 = openmc.Cell(
        name="FCC-FO-3", fill=m_graphite, region=-s_fcc_g_3 & +s_fcc_f_3
    )
    fcc_fo_4 = openmc.Cell(
        name="FCC-FO-4", fill=m_graphite, region=-s_fcc_g_4 & +s_fcc_f_4
    )
    fcc_fo_5 = openmc.Cell(
        name="FCC-FO-5", fill=m_graphite, region=-s_fcc_g_5 & +s_fcc_f_5
    )
    fcc_fo_6 = openmc.Cell(
        name="FCC-FO-6", fill=m_graphite, region=-s_fcc_g_6 & +s_fcc_f_6
    )
    fcc_fo_7 = openmc.Cell(
        name="FCC-FO-7", fill=m_graphite, region=-s_fcc_g_7 & +s_fcc_f_7
    )
    fcc_fo_8 = openmc.Cell(
        name="FCC-FO-8", fill=m_graphite, region=-s_fcc_g_8 & +s_fcc_f_8
    )
    fcc_fo_9 = openmc.Cell(
        name="FCC-FO-9", fill=m_graphite, region=-s_fcc_g_9 & +s_fcc_f_9
    )
    fcc_fo_10 = openmc.Cell(
        name="FCC-FO-10", fill=m_graphite, region=-s_fcc_g_10 & +s_fcc_f_10
    )
    fcc_fo_11 = openmc.Cell(
        name="FCC-FO-11", fill=m_graphite, region=-s_fcc_g_11 & +s_fcc_f_11
    )
    fcc_fo_12 = openmc.Cell(
        name="FCC-FO-12", fill=m_graphite, region=-s_fcc_g_12 & +s_fcc_f_12
    )
    fcc_fo_13 = openmc.Cell(
        name="FCC-FO-13", fill=m_graphite, region=-s_fcc_g_13 & +s_fcc_f_13
    )
    fcc_fo_14 = openmc.Cell(
        name="FCC-FO-14", fill=m_graphite, region=-s_fcc_g_14 & +s_fcc_f_14
    )

    fcc_he_inter2 = openmc.Cell(
        name="Helium-Inter",
        fill=all_coolant,
        region=+s_fcc_g_1
        & +s_fcc_g_2
        & +s_fcc_g_3
        & +s_fcc_g_4
        & +s_fcc_g_5
        & +s_fcc_g_6
        & +s_fcc_g_7
        & +s_fcc_g_8
        & +s_fcc_g_9
        & +s_fcc_g_10
        & +s_fcc_g_11
        & +s_fcc_g_12
        & +s_fcc_g_13
        & +s_fcc_g_14,
    )

    fcc_f_u_high_ = openmc.Universe(
        name="FCC-F-Raw",
        cells=[
            fcc_fi_1_high,
            fcc_fi_2_high,
            fcc_fi_3_high,
            fcc_fi_4_high,
            fcc_fi_5_high,
            fcc_fi_6_high,
            fcc_fi_7_high,
            fcc_fi_8_high,
            fcc_fi_9_high,
            fcc_fi_10_high,
            fcc_fi_11_high,
            fcc_fi_12_high,
            fcc_fi_13_high,
            fcc_fi_14_high,
            fcc_fo_1,
            fcc_fo_2,
            fcc_fo_3,
            fcc_fo_4,
            fcc_fo_5,
            fcc_fo_6,
            fcc_fo_7,
            fcc_fo_8,
            fcc_fo_9,
            fcc_fo_10,
            fcc_fo_11,
            fcc_fo_12,
            fcc_fo_13,
            fcc_fo_14,
            fcc_he_inter2,
        ],
    )
    fcc_f_c_high = openmc.Cell(
        name="FCC-Fuel-Pebbles", fill=fcc_f_u_high_, region=box_fcc_pebbles
    )
    fcc_f_u_high = openmc.Universe(name="FCC-Fuel-Pebbles", cells=[fcc_f_c_high])

    n_fuel = 8
    fuel_pebbles_8x8x10_lattice = openmc.RectLattice(name="FCC-8x8x10-Fuel-Pebbles")
    fuel_pebbles_8x8x10_lattice.lower_left = (
        -pebbles_pitch * n_fuel / 2,
        -pebbles_pitch * n_fuel / 2,
        -pebbles_pitch * 10 / 2,
    )
    fuel_pebbles_8x8x10_lattice.pitch = (pebbles_pitch, pebbles_pitch, pebbles_pitch)
    fuel_pebbles_8x8x10_lattice.outer = all_coolant
    fuel_pebbles_8x8x10_lattice.universes = [[[fcc_f_u_high] * n_fuel] * n_fuel] * 10

    pebble_bed_8x8x10_cell = openmc.Cell(
        name="PebbleBed-8x8x10-Cell", fill=fuel_pebbles_8x8x10_lattice
    )
    pebble_bed_8x8x10_universe = openmc.Universe(
        name="PebbleBed-8x8x10-Universe", cells=[pebble_bed_8x8x10_cell]
    )

    fuel_pebbles_lattice = openmc.RectLattice(name="FCC-Fuel-Pebbles")
    fuel_pebbles_lattice.lower_left = (-160, -160.0, 325 - lower_the_reactor)
    fuel_pebbles_lattice.pitch = (
        pebbles_pitch * n_fuel,
        pebbles_pitch * n_fuel,
        pebbles_pitch * 10,
    )
    fuel_pebbles_lattice.outer = all_coolant
    fuel_pebbles_lattice.universes = [
        [[pebble_bed_8x8x10_universe] * (int(40 / n_fuel) + 1)] * (int(40 / n_fuel) + 1)
    ] * 13

    return fuel_pebbles_lattice


print(f"Generation took {datetime.timedelta(seconds=time.perf_counter()-gen_time)}")


###============================================== ###
###==================== Refuel ================== ###
###============================================== ###
def get_coords_mat(index: int, model: openmc.Model) -> dict:
    results = openmc.deplete.Results(
        f"{cwd}/step_{index-1}/depletion_result_{index-1}.h5"
    )
    last_result = results[-1]
    geometry = model.geometry
    l_mat = geometry.get_lattices_by_name("FCC-Fuel-Pebbles")
    c_reg = geometry.get_cells_by_name("Pebble-Bed-Cell")
    c_reg = [
        f'{c.name.split("Pebble-Bed-Cell-")[1]}-{c.fill.id}'
        for c in c_reg
        if "FCC-Fuel-Pebbles" in c.fill.name
    ]
    coords_mat = {}
    for c in c_reg:
        coords, l_id = c.split("-")
        coords = int(coords)
        l_id = int(l_id)
        for l in l_mat:
            if l.id == l_id:
                u_tmp = openmc.Universe(cells=list(l.universes[0][0][0].cells.values()))
                g_tmp = openmc.Geometry(u_tmp)
                c_tmp = g_tmp.get_cells_by_name("Kernel")[0]
                m_tmp = c_tmp.fill
                mat = last_result.get_material(str(m_tmp.id))
                if fuel == "UCO":
                    mat.add_element("C", 9.92e-3, "ao")
                mat.temperature = coretemp
                coords_mat[coords] = [mat, int(c_tmp.name.split("Kernel-High-")[1])]
    return coords_mat


def refuel(coords_mat: dict):
    pb = []
    lowest_fuel = min(coords_mat.keys())
    highest_fuel = max(coords_mat.keys())
    for i in range(1, region + 1):
        if (
            highest_fuel < max_fuel
        ):  # Kalau belum, fuel ditambahkan di atas saja, tidak dikurangi
            if i < lowest_fuel:
                reg = openmc.Cell(
                    name=f"Pebble-Bed-Cell-{i}",
                    fill=graphite_pebbles_lattice,
                    region=+s_reg[i - 1] & -s_reg[i] & -s_r_2,
                )
            elif i >= lowest_fuel and i <= highest_fuel:
                fuel_lattice = build_fuel_lattice(coords_mat[i][0], coords_mat[i][1])
                reg = openmc.Cell(
                    name=f"Pebble-Bed-Cell-{i}",
                    fill=fuel_lattice,
                    region=+s_reg[i - 1] & -s_reg[i] & -s_r_2,
                )
            elif i == highest_fuel + 1:
                fuel_lattice = build_fuel_lattice(m_fuel_high)
                reg = openmc.Cell(
                    name=f"Pebble-Bed-Cell-{i}",
                    fill=fuel_lattice,
                    region=+s_reg[i - 1] & -s_reg[i] & -s_r_2,
                )
            else:
                reg = openmc.Cell(
                    name=f"Pebble-Bed-Cell-{i}",
                    fill=m_pendingin,
                    region=+s_reg[i - 1] & -s_reg[i] & -s_r_2,
                )
        else:  # Kalau fuel di core sudah banyak, baru refuel dengan mengganti fuel
            if i < lowest_fuel - 1:
                reg = openmc.Cell(
                    name=f"Pebble-Bed-Cell-{i}",
                    fill=graphite_pebbles_lattice,
                    region=+s_reg[i - 1] & -s_reg[i] & -s_r_2,
                )
            elif (
                i >= lowest_fuel - 1 and i <= highest_fuel - 1
            ):  # Lower Fuel by 1 Region
                fuel_lattice = build_fuel_lattice(
                    coords_mat[i + 1][0], coords_mat[i + 1][1]
                )
                reg = openmc.Cell(
                    name=f"Pebble-Bed-Cell-{i}",
                    fill=fuel_lattice,
                    region=+s_reg[i - 1] & -s_reg[i] & -s_r_2,
                )
            elif i == highest_fuel:
                if (
                    len(coords_mat.keys()) < max_fuel
                ):  # Add Fresh Fuel jika fuel belom penuh
                    mat = [m_fuel_high, 1]
                elif (
                    coords_mat[1][1] < regression
                ):  # shuffle if lowest fuel isnt at regression cap
                    mat = [coords_mat[1][0], coords_mat[1][1] + 1]
                else:  # Add fresh fuel jika regression sudah full
                    mat = [m_fuel_high, 1]
                fuel_lattice = build_fuel_lattice(mat[0], mat[1])
                reg = openmc.Cell(
                    name=f"Pebble-Bed-Cell-{i}",
                    fill=fuel_lattice,
                    region=+s_reg[i - 1] & -s_reg[i] & -s_r_2,
                )
            else:
                reg = openmc.Cell(
                    name=f"Pebble-Bed-Cell-{i}",
                    fill=m_pendingin,
                    region=+s_reg[i - 1] & -s_reg[i] & -s_r_2,
                )

        pb.append(reg)

    return pb


def build_geom(index: int, model: openmc.Model, pb: list):
    geometry = model.geometry
    old_pb = geometry.get_cells_by_name("Pebble-Bed-Cell")
    if len(pb) != len(old_pb):
        raise ValueError("Jumlah Pebble Bed Cell tidak Sesuai")
    root_universe = geometry.root_universe
    for i in range(len(old_pb)):
        root_universe.remove_cell(old_pb[i])
        root_universe.add_cell(pb[i])
    os.chdir(f"{cwd}/step_{index}")
    if test:
        new_tallies = None
    else:
        new_tallies = model.tallies
    new_geometry = openmc.Geometry(root_universe)
    new_materials = openmc.Materials(list(new_geometry.get_all_materials().values()))
    new_model = openmc.Model(new_geometry, new_materials, settings, new_tallies)
    os.chdir("../../")
    return new_model


###============================================== ###
###=================== Settings ================= ###
###============================================== ###
settings = openmc.Settings()
if test:
    settings.particles = 500
    settings.inactive = 25
    settings.batches = 100
else:
    settings.particles = 3000
    settings.inactive = 50
    settings.batches = 250

s_r_1 = openmc.ZCylinder(r=25.0)
s_r_2 = openmc.ZCylinder(r=25.0 + 125.275)
core_region = +s_reg[0] & -s_reg[region] & -s_r_1
lower_left, upper_right = core_region.bounding_box
uniform_dist = openmc.stats.Box(lower_left, upper_right)
source = openmc.IndependentSource(space=uniform_dist)
settings.source = source
settings.max_lost_particles = 100000

settings.temperature = {"method": "interpolation"}
settings.output = {"tallies": False}


###============================================== ###
###=================== Deplete ================== ###
###============================================== ###
def get_timesteps(X: int, N: int) -> list:
    base, rem = X // N, X % N
    return [base + 1] * rem + [base] * (N - rem)


# timesteps: list = get_timesteps(timestep, numberstep)
timesteps = [timestep / regression]
power = 250e6

chain = openmc.deplete.Chain.from_xml("../OpenMC/endfb71/chain_endfb71_pwr.xml")

os.makedirs(f"{cwd}/step_{startstep}", exist_ok=True)
print("Burnup begins")
cal_time = time.perf_counter()

if startstep != 0:
    coords_mat = get_coords_mat(startstep, model)
    pb = refuel(coords_mat)
    model = build_geom(startstep, model_new or model, pb)

operator = openmc.deplete.CoupledOperator(
    model, "../OpenMC/endfb71/chain_endfb71_pwr.xml"
)
integrator = openmc.deplete.PredictorIntegrator(
    operator, timesteps, power, timestep_units="s"
)
os.chdir(f"{cwd}/step_{startstep}")
settings.export_to_xml()
integrator.integrate(path=f"depletion_result_{startstep}.h5")
os.chdir("../../")
print(f"Calc took {datetime.timedelta(seconds=time.perf_counter()-cal_time)}")
