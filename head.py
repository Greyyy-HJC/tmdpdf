# %%
'''
Some useful settings
'''

grey = "#808080" 
red = "#FF6F6F" 
peach = "#FF9E6F" 
orange = "#FFBC6F" 
sunkist = "#FFDF6F"
yellow = "#FFEE6F"
lime = "#CBF169"
green = "#5CD25C" 
turquoise = "#4AAB89"
blue = "#508EAD" 
grape = "#635BB1"
violet = "#7C5AB8" 
fuschia = "#C3559F"

color_ls = ['orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue', 'red','green','orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue', 'red','green']

fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
fig_size  = (fig_width, fig_width / gr)
plt_axes = [0.12,0.15,0.8,0.78] # left, bottom, width, height
gridspec_tmin = {'height_ratios': [3, 1, 1], 'left': 0.15, 'right': 0.95, 'bottom': 0.16, 'top': 0.95} #* for tmin stability plot

errorp = {"markersize": 5, "mfc": "none", "linestyle": "none"} # circle
errorb = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1} # circle
errorl = {"markersize": 5, "mfc": "none", "capsize": 3, "elinewidth": 1} # circle with line
fs_p = {"fontsize": 13} # font size of text, label, ticks
ls_p = {"labelsize": 13}

font = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 15}

small_font = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 13}

gev_fm = 0.1973269631 # 1 = 0.197 GeV . fm

#* for latex to use Times New Roman
from matplotlib import rcParams
config = {
    "font.family":'serif',
    "font.size": 12,
    "mathtext.fontset":'stix',
    "font.serif": ['Times New Roman'],
}
rcParams.update(config)




#* Qi An's settings for tmdpdf project
def color_ls_Qi_An(scheme):
    if(scheme==1):
        colors = ['orange','dodgerblue','blueviolet','deeppink','royalblue','rosybrown', \
            'fuchsia','red','green','cyan', 'orange','dodgerblue',\
            'blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia',\
            'royalblue','red','green']
    if(scheme==2):
        colors = ["#4197d8", "#f8c120", "#413496", "#495226", "#d60b6f", "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d"]
    if(scheme==3):
        colors = ["#0055AA", "#C40003", "#00C19B", "#EAC862", "#7FD2FF", "#007ED3", "#B2DF8A", "#FFACAA", "#FF9D1E", "#C3EF00", "#CAB2D6", "#894FC6"]
    if(scheme==4):
        colors = ["#00441B", "#46A040", "#00AF99", "#FFC179", "#98D9E9", "#F6313E", "#FFA300", "#C390D4", "#FF5A00", "#AFAFAF", "#7D7D7D", \
            "#4B4B4B", "#8F1336", "#0081C9", "#001588", "#490C65", "#BA7FD0"]
    return colors

fig_size_Qi_An = (5, 3.09)
legend_fontsize_Qi_An = 10
fontsize_Qi_An = 12

font_Qi_An = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 12}