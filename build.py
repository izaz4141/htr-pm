# @title HTR-PM

import openmc
import openmc.model
import openmc.mgxs

import numpy as np
import math
import time
import datetime

mode = "deplete"
fuel = "UCO"
enrichment = "high"
test = True
wholemattemp = 273 + 20  # K
coretemp = 273 + 20  # Kelvin
seed = 666666
# Ditemukan bahwa jika partikel bereaksi dengan lattice pada koordinat tinggi, maka partikel tersebut akan hilang
lower_the_reactor = 1200
region = 60
max_fuel = 56
pebbles_pitch = 9.05164493902406
fcc_pebbles_pitch = pebbles_pitch / 2
if enrichment == 'low':
    packing_fraction = 0.072710770688
    triso_pitch = 0.177659234
else:
    packing_fraction = 0.25
    triso_pitch = 0.109392023
fuel_volume = 1167.33*3.14159*((25+125.275)**2)*0.61*packing_fraction*(2.5**3)*(250**3)/(3**3)/(460**3)/region
if mode == "deplete":
    cr_height = 100
else:
    cr_height = 0

tinggi_mixed_pebbles = 275.0

def coolant_density(t):
  p_in_bar = 3e+6 * 1.0e-5;
  return 48.14 * p_in_bar / (t + 0.4446 * p_in_bar / math.pow(t, 0.2)) / 10;


if __name__ == "__main__":

    gen_time = time.perf_counter()

    np.random.seed(seed)

    m_fuel = openmc.Material(name="Fuel-UO2-4,2%")
    m_fuel.add_nuclide("U235", 0.037022077, "wo")
    m_fuel.add_nuclide("U238", 0.844455942, "wo")
    m_fuel.add_element("O", 0.118517981, "wo")
    m_fuel.add_element("B", 4e-6, "wo")
    m_fuel.add_s_alpha_beta('c_U_in_UO2')
    m_fuel.set_density("g/cm3", 10.4)
    m_fuel.temperature = coretemp

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

    w_B = 0.795e-6 / 4
    # Buffer
    m_Buffer = openmc.Material(name="Buffer")
    m_Buffer.set_density("g/cm3", 1.1)
    m_Buffer.add_element('C', 1-w_B, 'wo')
    m_Buffer.add_element('B', w_B, 'wo')
    m_Buffer.add_s_alpha_beta("c_Graphite")
    m_Buffer.temperature = coretemp

    m_PyC = openmc.Material(name="m_PyC")
    m_PyC.set_density("g/cm3", 1.9)
    m_PyC.add_element('C', 1-w_B, 'wo')
    m_PyC.add_element('B', w_B, 'wo')
    m_PyC.add_s_alpha_beta("c_Graphite")
    m_PyC.temperature = coretemp

    m_SiC = openmc.Material(name="m_SiC")
    w_Si = 28.0855/40.0965 * (1 - w_B)
    w_C  = 12.011/40.0965  * (1 - w_B)
    m_SiC.add_element("Si", w_Si, "wo")
    m_SiC.add_element("C",  w_C,  "wo")
    m_SiC.add_element("B",  w_B,  "wo")
    m_SiC.add_s_alpha_beta("c_Graphite")
    m_SiC.set_density("g/cm3", 3.18)
    m_SiC.temperature = coretemp

    m_graphite = openmc.Material(name="m_Graphite_Matrix")
    m_graphite.set_density("g/cm3", 1.74)
    m_graphite.add_element('C', 1-0.795E-6, 'wo')
    m_graphite.add_element('B', 0.795E-6, 'wo')
    m_graphite.add_s_alpha_beta("c_Graphite")
    m_graphite.temperature = coretemp

    m_graphite_pebble = openmc.Material(name="m_Graphite_Pebble")
    m_graphite_pebble.set_density("g/cm3", 1.74)
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
        m_pendingin.add_element('H', 3.41344e-6, 'ao')
        m_pendingin.add_element('N', 3.44020e-5, 'ao')
        m_pendingin.add_element('O', 1.09688e-5, 'ao')
        m_pendingin.set_density('g/cm3', 1.092e-3)
    m_pendingin.temperature = coretemp

    # --------------------------------------------------------------------------------------------------

    m_cr = openmc.Material(name="m_Control-Rod")
    m_cr.add_element("C", 1.04589e-2, "ao")
    m_cr.add_nuclide("B10", 3.35103e-2, "ao")
    m_cr.add_nuclide("B11", 8.32529e-3, "ao")
    m_cr.add_element("Fe", 8.64966e-3, "ao")
    m_cr.add_element("He", 1.99240e-4, "ao")
    m_cr.temperature = wholemattemp
    m_cr.set_density("g/cm3", 1.7) # HTR-10 Value

    m_nbcb = openmc.Material(name="Non-Borated-CarbonBrick")
    m_nbcb.add_element("C", 8.53015e-2, "ao")
    m_nbcb.add_element("B", 3.7919e-5, "ao")
    m_nbcb.temperature = wholemattemp
    m_nbcb.add_s_alpha_beta("c_Graphite")
    m_nbcb.set_density("sum")

    m_bcb = openmc.Material(name="Borated-CarbonBrick")
    m_bcb.add_element("C", 8.39217e-2, "ao")
    m_bcb.add_element("B", 3.79675e-3, "ao")
    m_bcb.temperature = wholemattemp
    m_bcb.add_s_alpha_beta("c_Graphite")
    m_bcb.set_density("sum")

    IG110 = [8.92947e-2, 4.41429e-8, 1.781]
    m_sr = openmc.Material(name="Standard-Reflector")
    m_sr.add_element("C", IG110[0], "ao")
    m_sr.add_element("B", IG110[1], "ao")
    m_sr.temperature = wholemattemp
    m_sr.add_s_alpha_beta("c_Graphite")
    m_sr.set_density("g/cm3", IG110[2])

    m_reflector1 = openmc.Material(name="Reflector-52")
    m_reflector1.add_element("C", 0.5079 * IG110[0], "ao")
    # m_reflector1.add_element("B", 0.5079 * IG110[1], "ao")
    m_reflector1.temperature = wholemattemp
    m_reflector1.add_s_alpha_beta("c_Graphite")
    m_reflector1.set_density("g/cm3", IG110[2])

    m_reflector2 = openmc.Material(name="Reflector-58")
    m_reflector2.add_element("C", 0.9317 * IG110[0], "ao")
    # m_reflector2.add_element("B", 0.9317 * IG110[1], "ao")
    m_reflector2.temperature = wholemattemp
    m_reflector2.add_s_alpha_beta("c_Graphite")
    m_reflector2.set_density("g/cm3", IG110[2])

    m_reflector3 = openmc.Material(name="Reflector-55")
    m_reflector3.add_element("C", 0.8741 * IG110[0], "ao")
    # m_reflector3.add_element("B", 0.8741 * IG110[1], "ao")
    m_reflector3.temperature = wholemattemp
    m_reflector3.add_s_alpha_beta("c_Graphite")
    m_reflector3.set_density("g/cm3", IG110[2])

    m_reflector10 = openmc.Material(name="Reflector-46")
    m_reflector10.add_element("C", 0.6679 * IG110[0], "ao")
    # m_reflector10.add_element("B", 0.6679 * IG110[1], "ao")
    m_reflector10.temperature = wholemattemp
    m_reflector10.add_s_alpha_beta("c_Graphite")
    m_reflector10.set_density("g/cm3", IG110[2])

    m_brp1 = openmc.Material(name="Reflector-Poison-1")
    m_brp1.add_element("C", 7.49075e-2, "ao")
    m_brp1.add_element("B", 1.26898e-4, "ao")
    m_brp1.temperature = wholemattemp
    m_brp1.add_s_alpha_beta("c_Graphite")
    m_brp1.set_density("sum")

    m_brp2 = openmc.Material(name="Reflector-Poison-2")
    m_brp2.add_element("C", 8.14783e-2, "ao")
    m_brp2.add_element("B", 1.22899e-4, "ao")
    m_brp2.temperature = wholemattemp
    m_brp2.add_s_alpha_beta("c_Graphite")
    m_brp2.set_density("sum")

    m_reflector4 = openmc.Material(name="Reflector-49")
    m_reflector4.add_element("C", 0.7135 * IG110[0], "ao")
    # m_reflector4.add_element("B", 0.7135 * IG110[1], "ao")
    m_reflector4.temperature = wholemattemp
    m_reflector4.add_s_alpha_beta("c_Graphite")
    m_reflector4.set_density("g/cm3", IG110[2])

    m_reflector5 = openmc.Material(name="Reflector-54")
    m_reflector5.add_element("C", 0.8315 * IG110[0], "ao")
    # m_reflector5.add_element("B", 0.8315 * IG110[1], "ao")
    m_reflector5.temperature = wholemattemp
    m_reflector5.add_s_alpha_beta("c_Graphite")
    m_reflector5.set_density("g/cm3", IG110[2])

    m_reflector9 = openmc.Material(name="Reflector-19")
    m_reflector9.add_element("C", 0.719 * IG110[0], "ao")
    # m_reflector9.add_element("B", 0.719 * IG110[1], "ao")
    m_reflector9.temperature = wholemattemp
    m_reflector9.add_s_alpha_beta("c_Graphite")
    m_reflector9.set_density("g/cm3", IG110[2])

    m_reflector13 = openmc.Material(name="Reflector-38")
    m_reflector13.add_element("C", 0.99 * IG110[0], "ao")
    # m_reflector13.add_element("B", 0.99 * IG110[1], "ao")
    m_reflector13.temperature = wholemattemp
    m_reflector13.add_s_alpha_beta("c_Graphite")
    m_reflector13.set_density("g/cm3", IG110[2])

    m_reflector6 = openmc.Material(name="Reflector-3")
    m_reflector6.add_element("C", 0.91 * IG110[0], "ao")
    # m_reflector6.add_element("B", 0.91 * IG110[1], "ao")
    m_reflector6.temperature = wholemattemp
    m_reflector6.add_s_alpha_beta("c_Graphite")
    m_reflector6.set_density("g/cm3", IG110[2])

    m_reflector7 = openmc.Material(name="Reflector-5")
    m_reflector7.add_element("C", 0.84 * IG110[0], "ao")
    # m_reflector7.add_element("B", 0.84 * IG110[1], "ao")
    m_reflector7.temperature = wholemattemp
    m_reflector7.add_s_alpha_beta("c_Graphite")
    m_reflector7.set_density("g/cm3", IG110[2])

    m_reflector11 = openmc.Material(name="Reflector-9")
    m_reflector11.add_element("C", 5.62021e-2, "ao")
    m_reflector11.add_element("B", 3.69779e-6, "ao")
    m_reflector11.temperature = wholemattemp
    m_reflector11.add_s_alpha_beta("c_Graphite")
    m_reflector11.set_density("sum")

    m_reflector8 = openmc.Material(name="Reflector-6")
    m_reflector8.add_element("C", 0.6286 * IG110[0], "ao")
    # m_reflector8.add_element("B", 0.6286 * IG110[1], "ao")
    m_reflector8.temperature = wholemattemp
    m_reflector8.add_s_alpha_beta("c_Graphite")
    m_reflector8.set_density("g/cm3", IG110[2])

    m_reflector12 = openmc.Material(name="Reflector-10")
    m_reflector12.add_element("C", 2.30380e-2, "ao")
    m_reflector12.add_element("B", 3.68139e-6, "ao")
    m_reflector12.temperature = wholemattemp
    m_reflector12.add_s_alpha_beta("c_Graphite")
    m_reflector12.set_density("sum")

    mat_time = time.perf_counter()
    print(
        f"---------------------Material Generated in {datetime.timedelta(seconds=(mat_time-gen_time))}---------------------"
    )

    # --------------------------------------------------------------------------------------------------
    if enrichment != "high":
        # UO2 TRISO
        kernel_r = openmc.Sphere(r=250e-4)
        buffer_r = openmc.Sphere(r=(250 + 95) * 1e-4)
        IPyC_r = openmc.Sphere(r=(250 + 95 + 40) * 1e-4)
        SiC_r = openmc.Sphere(r=(250 + 95 + 40 + 35) * 1e-4)
        OPyC_r = openmc.Sphere(r=(250 + 95 + 40 + 35 + 40) * 1e-4)
        r_triso = (250 + 95 + 40 + 35 + 40) * 1e-4

    else:
        # UCO TRISO
        kernel_r = openmc.Sphere(r=212.5e-4)
        buffer_r = openmc.Sphere(r=(212.5 + 100) * 1e-4)
        IPyC_r = openmc.Sphere(r=(212.5 + 100 + 40) * 1e-4)
        SiC_r = openmc.Sphere(r=(212.5 + 100 + 40 + 35) * 1e-4)
        OPyC_r = openmc.Sphere(r=(212.5 + 100 + 40 + 35 + 40) * 1e-4)
        r_triso = (212.5 + 100 + 40 + 35 + 40) * 1e-4

    if mode == "deplete":
        c_kernel = openmc.Cell(name="Kernel-High-4", fill=m_fuel_high, region=-kernel_r)
    else:
        c_kernel = openmc.Cell(name="Kernel-High-4", fill=m_fuel, region=-kernel_r)

    triso_particle = openmc.Universe(
        name="Triso-Particle",
        cells=[
            c_kernel,
            openmc.Cell(name="Buffer", fill=m_Buffer, region=+kernel_r & -buffer_r),
            openmc.Cell(name="IPyC", fill=m_PyC, region=+buffer_r & -IPyC_r),
            openmc.Cell(name="m_SiC", fill=m_SiC, region=+IPyC_r & -SiC_r),
            openmc.Cell(name="OPyC", fill=m_PyC, region=+SiC_r & -OPyC_r),
            openmc.Cell(name="All-Graphite", fill=m_graphite, region=+OPyC_r),
        ],
    )

    # pf usually 0.05 this one is 0.0727108
    # jumlah triso di HTR-PM 11672
    # pf UCO adalah 40% atau 25% (AGR-5/6/7)
    fuel_pebble_r = openmc.Sphere(r=2.5)
    pebble_r = openmc.Sphere(
        r=3, 
        boundary_type="reflective"
    )
    fuel_pebble = openmc.Cell(name="Fuel-Pebble", region=-fuel_pebble_r)

    # Lattice Triso Random Makan Waktu Banyak
    # centers = openmc.model.pack_spheres(radius=(250+95+40+35+40)*1e-4, region=-fuel_pebble_r, num_spheres=11672, seed=12352653725)
    # trisos = [openmc.model.TRISO((250+95+40+35+40)*1e-4,triso_particle, center) for center in centers]

    # lower_left, upper_right = fuel_pebble.region.bounding_box
    # shape = (1, 1, 1)
    # pitch = (upper_right - lower_left)/shape
    # triso_lattice = openmc.model.create_triso_lattice(
    #     trisos, lower_left, pitch, shape, m_graphite)
    # ------------------------------------------

    graphite_matrix = openmc.Universe(
        name="Graphite-Matrix",
        cells=[openmc.Cell(name="Graphite-Matrix", fill=m_graphite)],
    )
    all_coolant = openmc.Universe(
        name="All-Coolant-Universe",
        cells=[openmc.Cell(name="All-Coolant-Cell", fill=m_pendingin)],
    )

    # Lattice Triso FCC (Satu Universe 4 Triso Particle) 4.02997 65.44937
    # triso_pitch = 0.282016
    # fcc_pitch = 0.282016/2
    # top_f = openmc.ZPlane(fcc_pitch)
    # bottom_f = openmc.ZPlane(-fcc_pitch)
    # right_f = openmc.XPlane(fcc_pitch)
    # left_f = openmc.XPlane(-fcc_pitch)
    # front_f = openmc.YPlane(fcc_pitch)
    # behind_f = openmc.YPlane(-fcc_pitch)
    # box_f = -top_f & +bottom_f & -right_f & +left_f & -front_f & +behind_f

    # fcc_lattice_triso = openmc.RectLattice(name='FCC-Triso')
    # fcc_lattice_triso.lower_left = (-3*fcc_pitch/2,-3*fcc_pitch/2,-3*fcc_pitch/2)
    # fcc_lattice_triso.pitch = (fcc_pitch,fcc_pitch,fcc_pitch)
    # fcc_lattice_triso.outer = openmc.Universe(cells=[openmc.Cell(fill=m_graphite)])
    # fcc_lattice_triso.universes = [
    #     [[triso_particle, graphite_matrix, triso_particle],[graphite_matrix, triso_particle, graphite_matrix],[triso_particle, graphite_matrix, triso_particle]],
    #     [[graphite_matrix, triso_particle, graphite_matrix],[triso_particle, graphite_matrix, triso_particle],[graphite_matrix, triso_particle, graphite_matrix]],
    #     [[triso_particle, graphite_matrix, triso_particle],[graphite_matrix, triso_particle, graphite_matrix],[triso_particle, graphite_matrix, triso_particle]]
    # ]
    # fcc_cell = openmc.Cell(name='FCC-Triso', fill=fcc_lattice_triso, region=box_f)
    # fcc_universe = openmc.Universe(name='FCC-Triso', cells=[fcc_cell])

    # triso_lattice = openmc.RectLattice(name='FCC-Pebble')
    # triso_lattice.lower_left = (-3,-3,-3)
    # triso_lattice.pitch = (triso_pitch,triso_pitch,triso_pitch)
    # triso_lattice.outer = openmc.Universe(cells=[openmc.Cell(fill=m_graphite)])
    # triso_lattice.universes = [
    #     [ [ fcc_universe for i in range(24) ] for i in range(24) ] for i in range(24)
    # ]
    def circle_points(r, n):
        theta = np.linspace(0, 2*np.pi, n, endpoint=False) + np.pi/(n/2)
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        return np.column_stack((x, y))

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
    # --------------------------------------------------------------------------------------------------

    # Lattice Triso Satu kubus satu triso
    

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

    # -----------------------------------------------------------------------------------------------------

    fuel_pebble.fill = triso_lattice

    outer_graphite_pebble = openmc.Cell(
        name="Graphite-Outer-Pebble", fill=m_graphite, region=+fuel_pebble_r & -pebble_r
    )
    luar_pebble1 = openmc.Cell(name="Moderator-1", fill=m_pendingin, region=+pebble_r)
    fuel_pebble_universe = openmc.Universe(
        name="Fuel-Pebble", cells=[fuel_pebble, outer_graphite_pebble, luar_pebble1]
    )

    # graphite_pebble = openmc.Cell(name='Graphite-Pebble', fill=m_graphite, region=-pebble_r)
    # luar_pebble2 = openmc.Cell(name='Moderator-1', fill=m_pendingin, region=+pebble_r)
    # graphite_pebble_universe = openmc.Universe(name='Graphite-Pebble', cells=[graphite_pebble, luar_pebble2])

    # -----------------------------------------------------------------------------------------------------
    # Z Axis
    s_z_0 = openmc.ZPlane(0.0 - lower_the_reactor, boundary_type="vacuum")
    s_z_1 = openmc.ZPlane(40.0 - lower_the_reactor)
    s_z_2 = openmc.ZPlane(80.0 - lower_the_reactor)
    s_z_3 = openmc.ZPlane(100.0 - lower_the_reactor)
    s_z_4 = openmc.ZPlane(180.0 - lower_the_reactor)
    s_z_5 = openmc.ZPlane(200.0 - lower_the_reactor)
    s_z_6 = openmc.ZPlane(220.0 - lower_the_reactor)
    s_z_7 = openmc.ZPlane(326.17 - lower_the_reactor)

    # Initial Graphite Pebble Height
    s_z_8 = openmc.ZPlane(931.17 - lower_the_reactor)
    tinggi_upper_cavity = 562.33 - tinggi_mixed_pebbles

    h_region = (605 + 562.33) / region
    s_reg = [
        openmc.ZPlane(326.17 - lower_the_reactor + (i * h_region))
        for i in range(0, region + 1)
    ]
    s_z_9 = openmc.ZPlane(931.17 + tinggi_mixed_pebbles - lower_the_reactor)

    s_z_10 = openmc.ZPlane(1493.5 - lower_the_reactor)
    s_z_11 = openmc.ZPlane(1542.5 - lower_the_reactor)
    s_z_12 = openmc.ZPlane(1577.5 - lower_the_reactor)
    s_z_13 = openmc.ZPlane(1630 - lower_the_reactor)
    s_z_14 = openmc.ZPlane(1670.0 - lower_the_reactor, boundary_type="vacuum")

    # XY Axis
    s_r_1 = openmc.ZCylinder(r=25.0)
    s_r_2 = openmc.ZCylinder(r=25.0 + 125.275)
    s_r_3 = openmc.ZCylinder(r=25.0 + 125.275 + 5.625)
    s_r_4 = openmc.ZCylinder(r=25.0 + 125.275 + 5.625 + 13.2)
    s_r_5 = openmc.ZCylinder(r=25.0 + 125.275 + 5.625 + 13.2 + 7.8)
    s_r_6 = openmc.ZCylinder(r=25.0 + 125.275 + 5.625 + 13.2 + 7.8 + 7.6)
    s_r_7 = openmc.ZCylinder(r=25.0 + 125.275 + 5.625 + 13.2 + 7.8 + 7.6 + 7.4)
    s_r_8 = openmc.ZCylinder(r=25.0 + 125.275 + 5.625 + 13.2 + 7.8 + 7.6 + 7.4 + 18.2)
    s_r_9 = openmc.ZCylinder(
        r=25.0 + 125.275 + 5.625 + 13.2 + 7.8 + 7.6 + 7.4 + 18.2 + 15.313
    )
    s_r_10 = openmc.ZCylinder(
        r=25.0 + 125.275 + 5.625 + 13.2 + 7.8 + 7.6 + 7.4 + 18.2 + 15.313 + 25.046,
        boundary_type="vacuum",
    )  # 250,459 cm

    # Control Rod Channel
    s_z_cr = openmc.ZPlane(1493.5 - lower_the_reactor - cr_height)
    coords_cr = circle_points(25.0 + 125.275 + 5.625 + 6.6, 30)
    coords_cr = np.delete(coords_cr, np.arange(0, 30, 5), axis=0)
    s_cr = [openmc.ZCylinder(r=6.5, x0=x, y0=y) for x, y in coords_cr]
    # Cold Helium Channel
    coords_ch = circle_points(25.0 + 125.275 + 5.625 + 13.2 + 7.8 + 7.6 + 7.4 + 9.1, 30)
    s_ch = [openmc.ZCylinder(r=9, x0=x, y0=y) for x, y in coords_ch]

    # -----------------------------------------------------------------------------------------------------

    # FCC GRAPHITE PEBBLE

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

    # FCC FUEL PEBBLE
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

    # Inner Part of Pebble (Fuel)
    fcc_fi_1 = openmc.Cell(name="FCC-FI-1", fill=triso_lattice, region=-s_fcc_f_1)
    fcc_fi_2 = openmc.Cell(name="FCC-FI-2", fill=triso_lattice, region=-s_fcc_f_2)
    fcc_fi_3 = openmc.Cell(name="FCC-FI-3", fill=triso_lattice, region=-s_fcc_f_3)
    fcc_fi_4 = openmc.Cell(name="FCC-FI-4", fill=triso_lattice, region=-s_fcc_f_4)
    fcc_fi_5 = openmc.Cell(name="FCC-FI-5", fill=triso_lattice, region=-s_fcc_f_5)
    fcc_fi_6 = openmc.Cell(name="FCC-FI-6", fill=triso_lattice, region=-s_fcc_f_6)
    fcc_fi_7 = openmc.Cell(name="FCC-FI-7", fill=triso_lattice, region=-s_fcc_f_7)
    fcc_fi_8 = openmc.Cell(name="FCC-FI-8", fill=triso_lattice, region=-s_fcc_f_8)
    fcc_fi_9 = openmc.Cell(name="FCC-FI-9", fill=triso_lattice, region=-s_fcc_f_9)
    fcc_fi_10 = openmc.Cell(name="FCC-FI-10", fill=triso_lattice, region=-s_fcc_f_10)
    fcc_fi_11 = openmc.Cell(name="FCC-FI-11", fill=triso_lattice, region=-s_fcc_f_11)
    fcc_fi_12 = openmc.Cell(name="FCC-FI-12", fill=triso_lattice, region=-s_fcc_f_12)
    fcc_fi_13 = openmc.Cell(name="FCC-FI-13", fill=triso_lattice, region=-s_fcc_f_13)
    fcc_fi_14 = openmc.Cell(name="FCC-FI-14", fill=triso_lattice, region=-s_fcc_f_14)

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

    fcc_f_u_ = openmc.Universe(
        name="FCC-F-Raw",
        cells=[
            fcc_fi_1,
            fcc_fi_2,
            fcc_fi_3,
            fcc_fi_4,
            fcc_fi_5,
            fcc_fi_6,
            fcc_fi_7,
            fcc_fi_8,
            fcc_fi_9,
            fcc_fi_10,
            fcc_fi_11,
            fcc_fi_12,
            fcc_fi_13,
            fcc_fi_14,
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
    fcc_f_c = openmc.Cell(
        name="FCC-Fuel-Pebbles", fill=fcc_f_u_, region=box_fcc_pebbles
    )
    fcc_f_u = openmc.Universe(name="FCC-Fuel-Pebbles", cells=[fcc_f_c])
    isi_mixed = np.array([fcc_f_u, fcc_g_u])

    # Refueling Hole
    # Isinya Graphite
    refuel_pebbles_lattice = openmc.RectLattice(name="FCC-Refuel-Pebbles")
    refuel_pebbles_lattice.lower_left = (-26, -26, -1 - lower_the_reactor)
    refuel_pebbles_lattice.pitch = (pebbles_pitch, pebbles_pitch, pebbles_pitch)
    refuel_pebbles_lattice.outer = all_coolant
    refuel_pebbles_lattice.universes = [[[fcc_g_u] * 6] * 6] * 40

    if mode == "deplete":
        nbcb_region = -s_r_10 & -s_z_1 & +s_z_0
        bcb_region = -s_r_10 & -s_z_2 & +s_z_1
        b_sr_1_region = -s_r_7 & +s_z_2 & -s_z_3
        brhhc_region = +s_z_3 & -s_z_4 & -s_r_3
        brhhcl1_region = +s_z_4 & -s_z_5 & -s_r_2
        brp1_region = +s_z_5 & -s_z_6 & -s_r_2
        brhhcl2_region = +s_z_6 & -s_z_7 & -s_r_2
    else:
        rh = openmc.Cell(
            name="Refueling-Hole",
            fill=refuel_pebbles_lattice,
            region=+s_z_0 & -s_z_7 & -s_r_1,
        )
        nbcb_region = -s_r_10 & -s_z_1 & +s_z_0 & +s_r_1
        bcb_region = -s_r_10 & -s_z_2 & +s_z_1 & +s_r_1
        b_sr_1_region = -s_r_7 & +s_z_2 & -s_z_3 & +s_r_1
        brhhc_region = +s_z_3 & -s_z_4 & -s_r_3 & +s_r_1
        brhhcl1_region = +s_z_4 & -s_z_5 & -s_r_2 & +s_r_1
        brp1_region = +s_z_5 & -s_z_6 & -s_r_2 & +s_r_1
        brhhcl2_region = +s_z_6 & -s_z_7 & -s_r_2 & +s_r_1

    # Bottom-1
    # Non-borated Carbon Bricks
    nbcb = openmc.Cell(name="NBCB", fill=m_nbcb, region=nbcb_region)

    # Bottom-2
    # Borated Carbon Bricks
    bcb = openmc.Cell(name="BCB", fill=m_bcb, region=bcb_region)

    # Bottom-3
    # Standard Reflector
    bottom_sr1 = openmc.Cell(name="Bottom-SR1", fill=m_sr, region=b_sr_1_region)
    bottom_sr2 = openmc.Cell(
        name="Bottom-SR2", fill=m_sr, region=-s_r_9 & +s_r_8 & +s_z_2 & -s_z_3
    )  # Dibolongi dengan

    # Bottom-4
    # Bottom reﬂector with hot helium chamber
    brhhc = openmc.Cell(name="BRHHC", fill=m_reflector1, region=brhhc_region)

    # Bottom reﬂector with hot helium guide tube
    brhhgt1 = openmc.Cell(
        name="BRHHGT-1", fill=m_reflector2, region=+s_z_3 & -s_z_4 & +s_r_3 & -s_r_7
    )
    brhhgt2 = openmc.Cell(
        name="BRHHGT-2", fill=m_reflector2, region=+s_z_3 & -s_z_4 & +s_r_8 & -s_r_9
    )

    # Bottom-5
    # Bottom reﬂector with hot helium channel
    brhhcl1 = openmc.Cell(name="BRHHCL-1", fill=m_reflector3, region=brhhcl1_region)

    bottom_sr3 = openmc.Cell(
        name="Bottom-SR3", fill=m_sr, region=+s_z_4 & -s_z_5 & +s_r_2 & -s_r_7
    )
    bottom_sr4 = openmc.Cell(
        name="Bottom-SR4", fill=m_sr, region=+s_z_4 & -s_z_5 & +s_r_8 & -s_r_9
    )

    # Bottom reﬂector with cold helium channel
    brchcl = openmc.Cell(
        name="BRCHCL", 
        fill=m_reflector10, 
        region=+s_z_2 & -s_z_5 & +s_r_7 & -s_r_8
        & openmc.Intersection((+cyl for cyl in s_ch))
    )

    # Bottom-6
    # Bottom reﬂector with B4C
    brp1 = openmc.Cell(name="Bottom-Reﬂector-B4C-1", fill=m_brp1, region=brp1_region)
    brp2 = openmc.Cell(
        name="Bottom-Reﬂector-B4C-2",
        fill=m_brp2,
        region=+s_z_5 & -s_z_6 & +s_r_2 & -s_r_9
        & openmc.Intersection((+cyl for cyl in s_ch))
    )

    # Bottom-7
    # Bottom reﬂector with hot helium channel
    brhhcl2 = openmc.Cell(name="BRHHCL-2", fill=m_reflector4, region=brhhcl2_region)

    # Bottom reﬂector with control rod channel
    brcrc = openmc.Cell(
        name="BRCRC", fill=m_reflector5, region=+s_z_6 & -s_z_7 & +s_r_3 & -s_r_4
    )

    ll_time = time.perf_counter()
    print(
        f"---------------------Lower Layer Generated in {datetime.timedelta(seconds=(ll_time-mat_time))}---------------------"
    )

    # Middle-Layer
    # Graphite Pebbles
    graphite_pebbles_lattice = openmc.RectLattice(name="FCC-Graphite-Pebbles")
    graphite_pebbles_lattice.lower_left = (-160, -160, 325 - lower_the_reactor)
    graphite_pebbles_lattice.pitch = (pebbles_pitch, pebbles_pitch, pebbles_pitch)
    graphite_pebbles_lattice.outer = all_coolant
    graphite_pebbles_lattice.universes = [[[fcc_g_u] * 40] * 40] * 70

    if mode != "deplete":
        gp = openmc.Cell(
            name="Graphite-Pebbles-Cell",
            fill=graphite_pebbles_lattice,
            region=+s_z_7 & -s_z_8 & -s_r_2,
        )

        # Mixed Pebbles & Upper Cavity
        mixed_pebbles_lattice = openmc.RectLattice(name="FCC-Mixed-Pebbles")
        mixed_pebbles_lattice.lower_left = (-160, -160.0, 325 - lower_the_reactor)
        mixed_pebbles_lattice.pitch = (pebbles_pitch, pebbles_pitch, pebbles_pitch)
        mixed_pebbles_lattice.outer = all_coolant
        mixed_pebbles_lattice.universes = [
            [
                [
                    np.random.choice(isi_mixed, 1, p=[0.46667, 0.53333])[0]
                    for x in range(50)
                ]
                for y in range(50)
            ]
            for z in range(170)
        ]
        mp = openmc.Cell(
            name="Mixed-Pebbles-Cell",
            fill=mixed_pebbles_lattice,
            region=+s_z_8 & -s_z_9 & -s_r_2,
        )
        uc = openmc.Cell(
            name="Upper-Cavity-Cell", fill=m_pendingin, region=+s_z_9 & -s_z_10 & -s_r_2
        )

    else:
        pb = []
        regression_pattern = [4,2,3,1]
        for i in range(1, region + 1):
            if i <= max_fuel:
                pattern_index = (i - 1) % len(regression_pattern)
                fuel_lattice = build_fuel_lattice(m_fuel_high.clone(), regression_pattern[pattern_index])
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

    # Inner Reflector
    ir = openmc.Cell(
        name="Inner-Reflector", fill=m_sr, region=+s_z_6 & -s_z_11 & +s_r_2 & -s_r_3
    )

    # Side reﬂector with control rod channel
    srcrc = openmc.Cell(
        name="SRCRC",
        fill=m_reflector9,
        region=+s_z_7 & -s_z_10 & +s_r_3 & -s_r_4 
        & openmc.Intersection((+cyl for cyl in s_cr))
    )

    # Middle Reflector
    mr1 = openmc.Cell(
        name="Middle-Reflector-1", fill=m_sr, region=+s_z_6 & -s_z_11 & +s_r_4 & -s_r_5
    )
    mr2 = openmc.Cell(
        name="Middle-Reflector-2", fill=m_sr, region=+s_z_6 & -s_z_11 & +s_r_6 & -s_r_7
    )

    # Side reﬂector with gap
    srg = openmc.Cell(
        name="SRG", fill=m_reflector13, region=+s_z_6 & -s_z_11 & +s_r_5 & -s_r_6
    )

    # Side reﬂector with cold helium channel
    srchcl = openmc.Cell(
        name="SRCHCL", 
        fill=m_reflector10, 
        region=+s_z_6 & -s_z_11 & +s_r_7 & -s_r_8
        & openmc.Intersection((+cyl for cyl in s_ch))
    )

    # Outer Reflector
    or1 = openmc.Cell(
        name="Outer-Reflector", fill=m_sr, region=+s_z_6 & -s_z_12 & +s_r_8 & -s_r_9
    )

    ml_time = time.perf_counter()
    print(
        f"---------------------Middle Layer Generated in {datetime.timedelta(seconds=(ml_time-mat_time))}---------------------"
    )

    # Top-Layer
    # Top reﬂector with charge tube
    trct = openmc.Cell(
        name="TRCT", fill=m_reflector6, region=+s_z_10 & -s_z_13 & -s_r_1
    )

    # Top reﬂector with cold helium channel
    trchcl = openmc.Cell(
        name="TRCHCL", fill=m_reflector7, region=+s_z_10 & -s_z_11 & +s_r_1 & -s_r_2
    )

    # Top reﬂector with control rod structure
    trcrs1 = openmc.Cell(
        name="TRCRS-1", 
        fill=m_reflector11, 
        region=+s_z_10 & -s_z_11 & +s_r_3 & -s_r_4
    )

    # Top reﬂector with cold helium chamber
    trchc1 = openmc.Cell(
        name="TRCHC1", fill=m_reflector8, region=+s_z_11 & -s_z_12 & +s_r_1 & -s_r_3
    ) 
    trchc2 = openmc.Cell(
        name="TRCHC2", 
        fill=m_reflector8, 
        region=+s_z_11 & -s_z_12 & +s_r_4 & -s_r_8
    ) 

    # Top reﬂector with control rod structure
    trcrs2 = openmc.Cell(
        name="TRCRS-2", 
        fill=m_reflector12, 
        region=+s_z_11 & -s_z_12 & +s_r_3 & -s_r_4
    )

    # Top Reflector
    tr1 = openmc.Cell(
        name="Top-Reflector-1", fill=m_sr, region=+s_z_12 & -s_z_13 & +s_r_1 & -s_r_3
    ) 
    tr2 = openmc.Cell(
        name="Top-Reflector-2", 
        fill=m_sr, 
        region=+s_z_12 & -s_z_13 & +s_r_4 & -s_r_9
    ) 

    # Top reﬂector with control rod structure
    trcrs3 = openmc.Cell(
        name="TRCRS-3", 
        fill=m_reflector11, 
        region=+s_z_12 & -s_z_13 & +s_r_3 & -s_r_4
    )

    # Control Rod
    control_rods = []
    empty_cr = []
    for i, cyl in enumerate(s_cr):
        control_rods.append(
            openmc.Cell(
                name=f"Control-Rod-{i}",
                fill=m_cr,
                region=+s_z_cr & -s_z_10 & -cyl
            )
        )
        empty_cr.append(
            openmc.Cell(
                name=f"Empty-Rod-{i}",
                fill=m_pendingin,
                region=+s_z_7 & -s_z_cr & -cyl
            )
        )
    # Empty Cold Helium Channel
    empty_ch = []
    for i, cyl in enumerate(s_ch):
        empty_ch.append(
            openmc.Cell(
                name=f"Empty-CHCL-{i}",
                fill=m_pendingin,
                region=+s_z_2 & -s_z_11 & -cyl
            )
        )


    # Outer Borated carbon bricks
    tobcb = openmc.Cell(
        name="Top-Borated-Carbon-Bricks", 
        fill=m_bcb, 
        region=+s_z_13 & -s_z_14 & -s_r_10
    )
    sobcb = openmc.Cell(
        name="Side-Borated-Carbon-Bricks",
        fill=m_bcb,
        region=+s_z_2 & -s_z_13 & +s_r_9 & -s_r_10,
    )

    tl_time = time.perf_counter()
    print(
        f"---------------------Top Layer Generated in {datetime.timedelta(seconds=(tl_time-mat_time))}---------------------"
    )
    univ_cells = [
        nbcb, bcb, bottom_sr1, bottom_sr2, brhhc, brhhgt1,
        brhhgt2, brhhcl1, bottom_sr3, bottom_sr4, brchcl, brp1,
        brp2, brhhcl2, brcrc, ir, srcrc, mr1,
        mr2, srg, srchcl, or1, trct, trchcl,
        trcrs1, trchc1, trchc2, trcrs2, tr1, tr2,
        trcrs3, tobcb, sobcb,
    ]
    univ_cells += control_rods
    univ_cells += empty_cr
    univ_cells += empty_ch
    if mode == "deplete":
        univ_cells += pb
    else:
        univ_cells += [rh, gp, mp, uc]

    main_universe = openmc.Universe(name="root_universe", cells=univ_cells)
    geometry = openmc.Geometry(main_universe)
    geometry.export_to_xml()

    materials = openmc.Materials(list(geometry.get_all_materials().values()))
    materials.export_to_xml()
    ex_time = time.perf_counter()
    print(
        f"--------------------- Exported in {datetime.timedelta(seconds=(ex_time - tl_time))}---------------------"
    )

    ### ============================================== ###
    ### =================== Settings ================= ###
    ### ============================================== ###
    settings = openmc.Settings()
    if test:
        settings.particles = 500
        settings.inactive = 50
        settings.batches = 200
    else:
        settings.particles = 10000
        settings.inactive = 50
        settings.batches = 200

    # lower_left = [-(25.+125.275+5.625+13.2+7.8+7.6+7.4+18.2+15.313+25.046), -(25.+125.275+5.625+13.2+7.8+7.6+7.4+18.2+15.313+25.046), 0]
    # upper_right = [(25.+125.275+5.625+13.2+7.8+7.6+7.4+18.2+15.313+25.046), (25.+125.275+5.625+13.2+7.8+7.6+7.4+18.2+15.313+25.046), (326.17+605+562.33+49.+35.+62.5+40.)]
    # lower_left, upper_right = fuel_pebble.region.bounding_box
    core_region = +s_reg[0] & -s_reg[region] & -s_r_1
    lower_left, upper_right = core_region.bounding_box
    uniform_dist = openmc.stats.Box(lower_left, upper_right)
    source = openmc.IndependentSource(space=uniform_dist)
    settings.source = source
    settings.max_lost_particles = 15000

    settings.temperature = {"method": "interpolation"}
    settings.output = {"tallies": False}
    settings.export_to_xml()

    ### ============================================== ###
    ### =================== Tallies ================== ###
    ### ============================================== ###
    # Instantiate an empty Tallies object
    tallies = openmc.Tallies()

    mesh = openmc.RegularMesh()
    mesh.dimension = [50, 50, 160]
    # mesh.lower_left = lower_left
    # mesh.upper_right = upper_right
    # mesh.lower_left, mesh.upper_right = mp.region.bounding_box # Dimensi Mixed Pebble
    mesh.lower_left, mesh.upper_right = main_universe.bounding_box

    normal_mesh = openmc.RegularMesh()
    normal_mesh.dimension = [1, 1, 1]
    normal_mesh.lower_left = main_universe.bounding_box[0]
    normal_mesh.upper_right = main_universe.bounding_box[1]

    # Create mesh filter for tally
    mesh_filter = openmc.MeshFilter(mesh)
    normal_filter = openmc.MeshFilter(normal_mesh)

    # Create energy filter for flux spectrum
    # energy_bins = [
    # 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 2.05, 5.0, 10.0, 50.0, 100.0, 500.0,
    # 1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7
    # ] # 20 Grup
    # SCALE -Group
    thermal = np.logspace(np.log10(1e-5), np.log10(0.625), num=6)
    epithermal = np.logspace(np.log10(0.625), np.log10(1e6), num=6)[1:]
    fast = np.logspace(np.log10(1e6), np.log10(20e8), num=3)[1:]
    energy_group = np.concatenate([thermal, epithermal, fast])

    energy_filter = openmc.EnergyFilter(energy_group)


    # Create Reaction Tally
    reaction_tally = openmc.Tally(name="Reaction")
    reaction_tally.filters = [mesh_filter]
    reaction_tally.nuclides = ["C0", "U235", "U238", "O16", "total"]
    reaction_tally.scores = [
        "fission",
        "(n,gamma)",
        "(n,a)",
        "absorption",
        "scatter",
        "total",
        "kappa-fission",
    ]

    # Create Spectrum Tally
    spectrum_tally = openmc.Tally(name="Spectrum")
    spectrum_tally.filters = [mesh_filter, energy_filter]
    spectrum_tally.scores = ["flux"]

    # Create Normalization Tally
    normal_tally = openmc.Tally(name="Heating-Local")
    normal_tally.filters = [normal_filter]
    normal_tally.scores = ["heating-local"]

    # Create beta tally
    beta = openmc.mgxs.Beta(
        name="Beta",
        domain=mesh,
        domain_type="mesh",
        energy_groups=openmc.mgxs.EnergyGroups(group_edges=energy_group),
        delayed_groups=list(range(1, 7)),
    )

    tallies.append(reaction_tally)
    tallies.append(spectrum_tally)
    tallies.append(normal_tally)
    tallies += beta.tallies.values()

    # Export to "tallies.xml"
    if not test:
        tallies.export_to_xml()

    print(f"Generation took {datetime.timedelta(seconds=time.perf_counter()-gen_time)}")
