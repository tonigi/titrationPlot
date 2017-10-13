
# coding: utf-8

# How to make SVG without X11: http://matplotlib.org/faq/howto_faq.html#matplotlib-in-a-web-application-server
# 
# How to embed in django: make SVG, replace in template http://stackoverflow.com/questions/34958702/embedding-a-matplotlib-plot-within-a-django-site
# 
# Alternatively, make PNG and embed as base64
# 

# In[2]:


get_ipython().magic('qtconsole')


# In[1]:


import pickle
import matplotlib
matplotlib.use('Agg')
import pandas as pd
from IPython.display import SVG, display


# In[2]:


pd=pickle.load(open("pd.pickle","rb"))


# In[4]:


# Alternatively
from htmd import *
pm,pd=proteinPrepare(Molecule("3PTB"),returnDetails=True)


# In[13]:


svg_plot=pd._get_pka_plot(pH=7.4)
print(svg_plot[0:400], "...")
display(SVG(svg_plot))


# In[24]:


def get_pka_plot(prepData, pH=7.4, figSizeX=10, dpk=1.0):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.colors import LinearSegmentedColormap
    
    
    # Shading
    Xe = np.array([[1,0],[1,0]])
    
    # Shading colors http://matplotlib.org/examples/pylab_examples/custom_cmap.html
    neutral_grey = (.7,.7,.7)
    my_red = (.98, .41, .29)
    my_blue = (.42, .68, .84)
    grey_red = LinearSegmentedColormap.from_list("grey_red",[neutral_grey,my_red ] )
    grey_blue = LinearSegmentedColormap.from_list("grey_blue",[neutral_grey,my_blue ] )
    eps = .01 # Tiny overprint to avoid very think white lines

    
    # Color for pk values
    pkcolor = "black"
    pkfontsize = 8
    dtxt = .3       # Displacement 
    
    
    # Or we could change the figure size, which scales axes
    # http://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
    SIZE = 12
    plt.rc('font', family="Open Sans")
    plt.rc('font', size=SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SIZE)  # legend fontsize
    plt.rc('figure', titlesize=SIZE)  # fontsize of the figure title

    
    # Constants
    acidicResidues = ['ASP','GLU','TYR']
    basicResidues = ['HIS','LYS','ARG']
    
    # titr =  (~ pd.isnull(d.pKa)) & d.pKa < 99
    d = prepData.data.copy()
    titr =  d.pKa < 99  # Automatically excludes NaN
    N = sum(titr)
    
    # Dubious residues
    d['dubious'] = abs(d.pKa - pH) < dpk
    
    # Format residue labels
    labels = ["{:s} {:s}:{:d}{:s}- {:s}".format("(!)" if x.dubious else "",
                                                x.chain,
                                                x.resid,
                                                x.insertion,
                                                x.resname) 
              for i,x in d.loc[titr].iterrows()]
    pKas = d.pKa.loc[titr]
    restypes = ["neg" if x.resname in acidicResidues else "pos" for i,x in d.loc[titr].iterrows() ]
    

    xmin, xmax = xlim = 0, 14
    ymin, ymax = ylim = -1, N
    
    width=.8    # Of each band
    
    # So, arbitrarily, 40 residues are square
    sizePerBand = figSizeX * (N/40)
    figsize=(figSizeX,sizePerBand)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, xlim=xlim, ylim=ylim,
                         autoscale_on=False)

    ax.xaxis.tick_top()
    ax.set_xlabel("pKa")
    ax.xaxis.set_label_position('top') 

    ax.yaxis.set_ticks(range(N))
    ax.yaxis.set_ticklabels(labels)
    ax.invert_yaxis()

        
    for i in range(N):
        left=xmin
        right=xmax
        top=i+width/2
        bottom=i-width/2
        pk = pKas.iloc[i]
        restype = restypes[i]
        
        if restype == "neg":
            ax.imshow(Xe*0, interpolation="none",
                      cmap=grey_blue, vmin=0, vmax=1,
                      extent=(left, pk-dpk, bottom, top), alpha=1)
            ax.imshow(np.fliplr(Xe), interpolation="bicubic",
                      cmap=grey_blue, vmin=0, vmax=1,
                      extent=(pk-dpk-eps, pk+dpk, bottom, top), alpha=1)
            ax.imshow(1+Xe*0, interpolation="none",
                      cmap=grey_blue, vmin=0, vmax=1,
                      extent=(pk+dpk-eps, right, bottom, top), alpha=1)
            ax.text(pk-dtxt,i," {:5.2f} ".format(pk),color=pkcolor, 
                    fontsize=pkfontsize, horizontalalignment="right",zorder=30)
        else:
            ax.imshow(1+Xe*0, interpolation="none",
                      cmap=grey_red, vmin=0, vmax=1,
                      extent=(left, pk-dpk, bottom, top), alpha=1)
            ax.imshow(Xe, interpolation="bicubic",
                      cmap=grey_red, vmin=0, vmax=1,
                      extent=(pk-dpk-eps, pk+dpk, bottom, top), alpha=1)
            ax.imshow(Xe*0, interpolation="none",
                      cmap=grey_red, vmin=0, vmax=1,
                      extent=(pk+dpk-eps, right, bottom, top), alpha=1)
            ax.text(pk+dtxt,i," {:5.2f} ".format(pk),color=pkcolor, 
                    fontsize=pkfontsize, horizontalalignment="left",zorder=30)
        ax.add_line(Line2D([pk,pk], [bottom,top], linewidth=3, color='white',zorder=2))

        # ax.add_line(Line2D([pk,pk], [bottom,top], linewidth=3, color='blue'))
        

    ## Shaded vertical band at pH  
    ax.axvline(x=pH-dpk, linewidth=2, color="black", alpha=.2, linestyle="dashed") 
    ax.axvline(x=pH+dpk, linewidth=2, color="black", alpha=.2, linestyle="dashed") 
    ax.axvline(x=pH, linewidth=3, color="black",alpha=.5) 
    ax.text(pH-dpk, ymax, " 90% protonated", rotation=90, 
            horizontalalignment="right", verticalalignment="bottom")
    ax.text(pH+dpk, ymax, " 10% protonated", rotation=90, 
            horizontalalignment="left", verticalalignment="bottom")
    

    ax.set_aspect('auto')
    
    # show()   # for interactive use
    from io import StringIO
    imgdata = StringIO()
    fig.savefig(imgdata, format="svg", bbox_inches='tight',)
    ret_img = imgdata.getvalue()
    
    fig.savefig("out.svg")
    fig.savefig("out.png") 
    
    # Png render may be a bit better -   
    # http://stackoverflow.com/questions/14824522/dynamically-serving-a-matplotlib-image-to-the-web-using-python
    ## 
    # from io import StringIO
    # buf = io.BytesIO()
    # plt.savefig(buf, format='png')
    # image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8').replace('\n', '')
    # buf.close()
    
    plt.close(fig)
    return ret_img


# In[25]:


svg_plot=get_pka_plot(pd, pH=7.4)
print(svg_plot[0:400], "...")
display(SVG(svg_plot))


# In[ ]:




