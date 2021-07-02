"""
pyscal visualization module
---------------------------

Used for visualization of system objects in jupyter notebooks and lab. Uses plotly.
"""
import numpy as np
from plotly import graph_objs as go
import ipywidgets as widgets
import itertools

def create_box_plot(box, origin=[0,0,0]):
    """
    Create traces which correspond to the simulation cell

    Parameters
    ----------
    box : list
        dimensions of the simulation box

    origin : list, optional
        Origin of the simulation box. Default [0, 0, 0]

    Returns
    -------
    traces : list of Scatter3d objects

    """
    box = np.array(box)
    origin = np.array(origin)
    combos = list(itertools.combinations(range(3), 2))
    faces = []
    for combo in combos:
        f1 = [origin, box[combo[0]], box[combo[0]]+box[combo[1]], box[combo[1]], origin]
        s = combo[0] + combo[1]
        t = 3-s
        f2 = [origin + box[t], box[combo[0]]+ box[t],  box[combo[0]]+box[combo[1]]+ box[t], box[combo[1]]+ box[t], origin + box[t]]
        faces.append(np.array(f1))
        faces.append(np.array(f2))
    traces = []
    for face in faces:
        trace = go.Scatter3d(
            x=face[:,0],
            y=face[:,1],
            z=face[:,2],
            mode='lines',
            name='lines',
            line=dict(width=2.0, color='#263238'),
            showlegend=False
        )
        traces.append(trace)
    return traces


def plot_3d(pos, color=None, radius=17, 
            colorscale='Spectral', opacity=1.0, 
            traces=None, cmap_title=None):
    """
    Plot the atoms along with the simulation box

    Parameters
    ----------
    pos : list of positions
        list of atomic positions

    color : list
        list of colors to use for plotting

    radius : int, optional
        radius of plotted atom objects

    colorscale : string, optional
        color map for coloring atoms

    opacity : float, optional
        opacity of atoms

    traces : box plot objects

    cmap_title : string
        title of cmap 
    """
    data=go.Scatter3d(
        x=pos[:,0],
        y=pos[:,1],
        z=pos[:,2],
        mode='markers',
        opacity=1.0,
        marker=dict(
            sizemode='diameter',
            sizeref=750,
            size=radius,
            color = color,
            opacity = opacity,
            colorscale = colorscale,
            colorbar=dict(thickness=20, title=cmap_title),
            line=dict(width=0.5, color='#455A64')
        )
    )
    traces.append(data)
    fig = go.Figure(data=traces)
    fig.update_layout(scene = dict(
                        xaxis_title="",
                        yaxis_title="",
                        zaxis_title="",
                        xaxis = dict(
                             showticklabels=False,
                             showbackground=False,
                             zerolinecolor="#455A64",),
                        yaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64"),
                        zaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64",),),
                        width=700,
                        margin=dict(
                        r=10, l=10,
                        b=10, t=10)
                      )
    fig.update_layout(showlegend=False)
    fig.show()

def plot_system(sys, colorby=None, filterby=None):
    """
    Plot the system

    Parameters
    ----------
    sys : System object

    colorby : string, optional
        property over which the atoms are to be colored. It can be any
        attributed of Atom, a custom attribute,  or calculated q values which can be accessed
        as `qx` or `aqx` where x stands for the q number.

    filterby : string, optional
        property over which the atoms are to be filtered before plotting.
        It can be any attribute of atom, or a custom value of atom. It should provide
        a True or False value.

    Returns
    -------
    None  
    """
    atoms = sys.atoms
    positions = []
    colors = []
    filters = []
    for count, atom in enumerate(atoms):
        if filterby is not None:
            if sys.get_custom(atom, [filterby])[0]:
                positions.append(np.array(sys.remap_atom(atom.pos)))
        else:
            positions.append(np.array(sys.remap_atom(atom.pos)))
    
        if colorby is not None:
            cx = sys.get_custom(atom, [colorby])[0]
            colors.append(cx)
        else:
            colors.append(1)
    colors = np.array(colors).astype(float)
    
    boxtraces = create_box_plot(sys.box)
    
    if colorby is None:
        ctitle = ""
    else:
        ctitle = colorby

    radius = widgets.FloatSlider(min=1, max=30, step=1)
    widgets.interact_manual.opts['manual_name'] = 'Render plot'
    im = widgets.interact_manual(plot_3d, radius=radius, 
        pos=widgets.fixed(np.array(positions)), cmap_title=widgets.fixed(ctitle),
        color=widgets.fixed(colors), opacity=widgets.fixed(1.0), 
        traces=widgets.fixed(boxtraces), description="test")