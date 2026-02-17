import openmc
import os, re, time, datetime
import numpy as np
import matplotlib.pyplot as plt

cwd = os.getcwd()
mode = 'normal'

def get_step_numbers():
    step_pattern = re.compile(r"^step_(\d+)$")  # Regex pattern to match "step_{number}"
    sim_pattern = "openmc_simulation_n"
    step_numbers = []
    sum_steps = 0
    
    for item in os.listdir(cwd):  # List items in the current directory with full path
        if os.path.isdir(os.path.join(cwd, item)):  # Check if it's a directory
            match = step_pattern.match(item)
            if match:
                step_numbers.append(int(match.group(1)))  # Extract number and convert to int
            for sim in os.listdir(f'{cwd}/{item}'):
                if sim_pattern in sim:
                    sum_steps += 1

    return sorted(step_numbers), sum_steps

def get_colors(geometry:openmc.Geometry):
    colors = {}
    graphites = geometry.get_materials_by_name('graphite')
    coolants = geometry.get_materials_by_name('pendingin')
    nbcb = geometry.get_materials_by_name('Non-Borated-CarbonBrick')[0]
    bcb = geometry.get_materials_by_name('Borated-CarbonBrick', matching=True)[0]
    poisones = geometry.get_materials_by_name('poison')
    control_rod = geometry.get_materials_by_name('reflector-19')[0]
    helium_channel = geometry.get_materials_by_name('reflector-46')[0]
    cr_top = geometry.get_materials_by_name('reflector-9')[0]
    cr_mid = geometry.get_materials_by_name('reflector-10')[0]
    sr = geometry.get_materials_by_name('standard-reflector')[0]
    colors[nbcb] = 'dimgray'
    colors[bcb] = 'olive'
    colors[control_rod] = 'lime'
    colors[helium_channel] = 'snow'
    colors[sr] = 'steelblue'
    colors[cr_top] = 'coral'
    colors[cr_mid] = 'coral'

    rids = [3, 5, 6, 49, 38, 54, 55, 52, 58]
    for rid in rids:
        reflector = geometry.get_materials_by_name(f'Reflector-{rid}', matching=True)[0]
        colors[reflector] = 'lightskyblue'

    for graphite in graphites:
        colors[graphite] = 'gray'
    for coolant in coolants:
        colors[coolant] = 'white'
    for poison in poisones:
        colors[poison] = 'green'
    return colors

step_numbers, sum_steps = get_step_numbers()

plot_time = time.perf_counter()
if mode == 'deplete':
    os.makedirs(f"{cwd}/geometry/", exist_ok=True)
    for step in step_numbers:
        model = openmc.Model.from_xml(
            f'{cwd}/step_{step}/geometry.xml',
            f'{cwd}/step_{step}/materials.xml',
            f'{cwd}/step_{step}/settings.xml',
            f'{cwd}/step_{step}/tallies.xml',
            f'{cwd}/step_{step}/plot.xml'
        )
        geometry = model.geometry
        colors = get_colors(geometry)
        

        print(f"##=================== Step: {step} ===================##")
        geometry.plot(color_by='material', 
                        pixels=(10000000),
                        # origin=(0,0,250),
                        colors=colors,
                        basis='xz')
        plt.savefig(f"{cwd}/geometry/reactor_xz_{step}.png")
        plt.close()

else:
    os.makedirs(f"{cwd}/eigenvalue/", exist_ok=True)
    model = openmc.Model.from_xml()
    geometry = model.geometry
    colors = get_colors(geometry)
    geometry.plot(color_by='cell', 
                    pixels=(100000000),
                    # width=(10,10),
                    origin=(0,0.05,835-1200),
                    # colors=colors,
                    basis='xz')
    plt.xlabel('x (cm)', fontsize=96)
    plt.ylabel('y (cm)', fontsize=96)
    plt.tick_params(axis='both', which='both', labelsize=96, length= 32, width= 7.5, direction= 'out')
    # plt.tight_layout()
    plt.savefig(f"{cwd}/eigenvalue/reactor_xz.png")
    plt.close()



# p.to_ipython_image()
# p = openmc.Plot.from_geometry(geometry)
# p.color_by = 'material'
# p.basis = 'xz'
# p.origin = (0,0,876 - lower_the_reactor)
# p.width = ((25.+125.275+5.625+13.2+7.8+7.6+7.4+18.2+15.313+25.046)*2,1670)
# p.colors = colors
# p.legend = True
# p.pixels = (int((25.+125.275+5.625+13.2+7.8+7.6+7.4+18.2+15.313+25.046)*2*5), int(1670*5))
# p.show_overlaps = True
# p.overlap_color =('snow',)
print(f"#================== Plot took {datetime.timedelta(seconds=time.perf_counter()-plot_time)} ==================#")