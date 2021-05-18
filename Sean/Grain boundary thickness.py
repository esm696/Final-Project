import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import plotly
import natsort
import plotly.graph_objects as go


def chunks(lst, n):
    """ breaks the list of data into chunks of fixed size """
    return [lst[x:x + n] for x in range(0, len(lst), n)]


def perimeter(dia):
    """perimeter of grains"""
    return np.pi * dia


# paths = []
# for path, dirs, files in os.walk(os.getcwd()):
# """another way to sort paths - set minimum length pf path as needed"""
#     #      print(len(path))
#     if len(path) >= 31:
#         paths.append(path)
# list_of_directories = natsort.os_sorted(paths)
# # print(paths)
# print(len(paths))

list_of_directories = []
# basic loop for directory generation
os.chdir('/Users/seanmascarenhas/Desktop/MoAu/MoAu2')
for i in range(50, 950, 50):
    for j in range(1, 31):
        list_of_directories.append(os.path.join(os.getcwd(), '{}/{}'.format(i, j)))
print(len(list_of_directories))

mean_size_list = []
d_list = []
std_list = []
fgb_list = []
t_sqrt = []
t_cbrt = []
p_list = []
sf_list=[]
unit = .314
height_sim = 4


def diameter(g_count):
    """diameter calculated from the number of atoms in the surface"""
    return (np.sqrt(4 * g_count/(height_sim * np.pi))) * unit

def shapefactor(area,peri):
    """shape factor for the grains"""
    return ((4*np.pi*area*unit*unit)/(height_sim*peri**2))


for items in list_of_directories:
    os.chdir(items)
    print(os.getcwd())
    step = pd.read_csv('step100000.txt', delim_whitespace='True', skiprows=9,
                       names=['x', 'y', 'z', 'cid', 'gid', 'gbid', 'nid'])

    grains = step.loc[step.gbid == 0]
    grains2 = grains[grains.duplicated(subset=['gid'], keep=False)]
    grains3 = grains2.groupby(['gid']).agg(counts_gid=('gid', 'count'))
    grains4 = grains3.loc[grains3.counts_gid >= 216].copy()
    grains4['diameter'] = grains4['counts_gid'].apply(lambda x: diameter(x))
    grains4['perimeter'] = grains4['diameter'].apply(lambda x: perimeter(x))
    grains4['shapefactor'] = grains4[['counts_gid','perimeter']].apply(lambda x: shapefactor(*x),axis=1)
    sf_list.append(grains4.shapefactor.mean())
    p = (600 * 600 * unit * unit) / grains4.perimeter.sum()
    p_list.append(p)
    mean_size_list.append(grains4['counts_gid'].mean())
    std_list.append(grains4['counts_gid'].std())
    avg_atoms_surf = grains4['counts_gid'].mean()
    d = diameter(avg_atoms_surf)
    d_list.append(d)
    fgb = step.gbid.sum()
    fgb = (fgb / 1440000)
    fgb_list.append(fgb)
    t_sqrt.append(d - (np.sqrt(1 - fgb)) * d)
    t_cbrt.append(d * (1 - (np.cbrt(1 - fgb))))

fig = go.FigureWidget(data=
go.Contour(
    z=chunks(t_sqrt, 30),
    y=np.array(range(50, 950, 50)),
    x=np.array(range(1, 31)),
    colorscale="viridis",
    colorbar=dict(
        title='GB thickness in nm',  # title here
        title_font_size=30,
        tickfont_size=20,
        titleside='right'
    )
    , contours=dict(
        start=0,
        end=6,
        size=.25,
        showlabels=True,  # show labels on contours
        labelfont=dict(  # label font properties
            size=18,
            color='white',
        )),
    #         contours_coloring='heatmap'
))
fig.update_layout(title_text='t  using sqrt shape factor and min size of 216',
                  xaxis_title='at. %',
                  yaxis_title='Temperature',
                  title=dict(font_size=20))
fig.update_xaxes(title=dict(font_size=30), tickfont_size=20)
fig.update_yaxes(title=dict(font_size=30), tickfont_size=20)
fig.update_layout(
    # autosize=False,
    # width=1200,
    # height=1000,
)

fig.show()

fig = go.FigureWidget(data=
go.Contour(
    z=chunks(t_cbrt, 30),
    y=np.array(range(50, 950, 50)),
    x=np.array(range(1, 31)),
    colorscale="ylgnbu",
    colorbar=dict(
        title='GB thickness in nm',  # title here
        title_font_size=30,
        tickfont_size=20,
        titleside='right'
    )
    , contours=dict(
        start=0,
        end=4.5,
        size=.25,
        showlabels=True,  # show labels on contours
        labelfont=dict(  # label font properties
            size=18,
            color='black',
        )),
    #         contours_coloring='heatmap'
))
fig.update_layout(title_text='t using cbrt and min size of 216',
                  xaxis_title='at. %',
                  yaxis_title='Temperature',
                  title=dict(font_size=20))
fig.update_xaxes(title=dict(font_size=30), tickfont_size=20)
fig.update_yaxes(title=dict(font_size=30), tickfont_size=20)

# fig.update_layout(
#     autosize=False,
#     width=1200,
#     height=1000,
# )

fig.show()

fig = go.FigureWidget(data=
go.Contour(
    z=chunks(p_list, 30),
    y=np.array(range(50, 950, 50)),
    x=np.array(range(1, 31)),
    colorscale="tempo",
    colorbar=dict(
        title='GB perimeter in nm',  # title here
        title_font_size=30,
        tickfont_size=20,
        titleside='right'
    ),
    contours=dict(
        # start=0,
        # end=6,
        # size=.25,
        showlabels=True,  # show labels on contours
        labelfont=dict(  # label font properties
            size=18,
            color='white',
        )),
    #         contours_coloring='heatmap'
))
fig.update_layout(title_text='Grain normalized perimeter with min cluster size of 216',
                  xaxis_title='at. %',
                  yaxis_title='Temperature',
                  title=dict(font_size=20))
fig.update_xaxes(title=dict(font_size=30), tickfont_size=20)
fig.update_yaxes(title=dict(font_size=30), tickfont_size=20)
fig.update_layout(
    # autosize=False,
    # width=1200,
    # height=1000,
)

fig.show()

fig = go.FigureWidget(data=
go.Contour(
    z=chunks(d_list, 30),
    y=np.array(range(50, 950, 50)),
    x=np.array(range(1, 31)),
    colorscale="portland",
    colorbar=dict(
        title='mean diameter in nm',  # title here
        title_font_size=30,
        tickfont_size=20,
        titleside='right'
    )
    , contours=dict(
        start=0,
        end=20,
        size=1,
        showlabels=True,  # show labels on contours
        labelfont=dict(  # label font properties
            size=18,
            color='white',
        )),
    #         contours_coloring='heatmap'
))
fig.update_layout(title_text='mean dia of grains using min cluster size of 216',
                  xaxis_title='at. %',
                  yaxis_title='Temperature',
                  title=dict(font_size=20))
fig.update_xaxes(title=dict(font_size=30), tickfont_size=20)
fig.update_yaxes(title=dict(font_size=30), tickfont_size=20)
# fig.update_layout(
#     autosize=False,
#     width=1200,
#     height=1000,
#
# )

fig.show()

fig = go.FigureWidget(data=
go.Contour(
    z=chunks(p_list, 30),
    y=np.array(range(50, 950, 50)),
    x=np.array(range(1, 31)),
    colorscale="solar",
    colorbar=dict(
        title='Effective length in nm',
        title_font_size=30,
        tickfont_size=20,
        titleside='right'
    )

    , contours=dict(
        start=-0,
        end=11.5,
        size=.25,
        showlabels=True,  # show labels on contours
        labelfont=dict(  # label font properties
            size=16,
            color='white',
        )),
    #             contours_coloring='heatmap'
))
fig.update_layout(title_text='Effective length and min size of 216',
                  xaxis_title='at. %',
                  yaxis_title='Temperature',
                  title=dict(font_size=20))
fig.update_xaxes(title=dict(font_size=30), tickfont_size=20)
fig.update_yaxes(title=dict(font_size=30), tickfont_size=20)

fig.update_layout(
    # autosize=False,
    # width=1200,
    # height=1000,
    #     margin=dict(
    #         l=50,
    #         r=50,
    #         b=100,
    #         t=100,
    #         pad=4
    #     ),
    #     paper_bgcolor="LightSteelBlue",
)

fig.show()