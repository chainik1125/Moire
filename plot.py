import numpy as np
import matplotlib.pyplot as plt
import json
import os
import ast
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.artist as artists
import matplotlib as mpl

from variables import *

labelfontsize=30
tickfontsize=30
ylabelsize=30
mpl.rcParams['text.usetex'] = True





def convert_to_json(data_dic):
    json_dic={}
    for i in data_dic.keys():
        json_dic[str(i)]=data_dic[i]
    return json_dic

#save data
def json_save(outputfilename,inputdic):
    with open(outputfilename, 'w') as fp:
        json.dump(inputdic, fp)

def json_open(filename):
    with open(filename, 'r') as fp:
        data = json.load(fp)
    return data

colors={1:'black',2:'darkorange',3:'gold',4:'darkgreen',5:'lime',6:'blue',7:'cornflowerblue',8:'red',9:'olive',10:'rebeccapurple',11:'darkslategrey',12:'maroon',13:'lightcoral',14:'chocolate',15:'navy',16:'magenta',17:'gold',18:'navy',19:'orange',20:'magenta',21:'darkblue',22:'maroon',23:'yellow',24:'red',25:'cyan',30:'brown',32:'maroon',33:'slateblue',34:'brown',35:'peru',36:'brown',48:'crimson',49:'purple',68:'yellow',70:'red',np.nan:'blue'}

directory = os.fsencode("/Users/dmitrymanning-coe/Documents/Research/Barry Bradlyn/Moire/Numerics/pathdata")
dstr="/Users/dmitrymanning-coe/Documents/Research/Barry Bradlyn/Moire/Numerics/pathdata"


def file_names(directory,contains,kpoints,variable,parameterv,theta):#omegas
    files=[]
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if (contains in filename) and (kpoints in filename) and (variable in filename) and (parameterv in filename) and (theta in filename):# 
            files.append(filename)
            continue
        else:
            continue
    if len(files)==0:
        print('No matches')
        return files
    else:
        print('More than one match')
        return files

# test_filenams=file_names(directory=directory,contains="AtoB",kpoints="5",variable="twoparticle",parameterv='UHK2')
# print(test_filenams)
# exit()

def group_parts(directory,dstr,start_string,end_string,kpoints,parts,variable,params,theta):
    files_list=file_names(directory=directory,dstr=dstr,contains=start_string+'to'+end_string,kpoints=kpoints,variable=variable,parameterv=params,theta=theta)
    
    values=[]
    for x in range(1,int(parts)+1):
        for xv in files_list:
            if f'part{x}of{parts}' in xv:
                temp=json_open(dstr+'/'+xv)
                if len(values)==0:
                    for rowindex in range(0,len(temp)):
                        rowtemp=temp[rowindex]
                        values.append(rowtemp)
                else:
                    for rowindex in range(0,len(temp)):
                        rowtemp=temp[rowindex]
                        for zv in rowtemp:
                            values[rowindex].append(zv)

    #You don't want to flatten completely. You want to add the rows in one
    return values

def group_parts_dic(directory,dstr,start_string,end_string,kpoints,params,variable,theta):
    files_list=file_names(directory=directory,contains=start_string+'to'+end_string,kpoints=kpoints,variable=variable,parameterv=params,theta=theta)
    #test_filenams=file_names(directory=directory2,contains="GammatoX",kpoints="500",variable="datadic")
    print(files_list)
    # print(files_list)
    # print(directory)
    # print(start_string)
    # print(end_string)
    # print(kpoints)

    # print(test_filenams)
    
    grouped_GF_dic={}
    # for x in range(1,int(parts)+1):
    #     for xv in files_list:
    #         if f'part{x}of{parts}' in xv:
    temp=json_open(dstr+'/'+files_list[0])
    for key in temp.keys():
        keyt=ast.literal_eval(key)
        grouped_GF_dic[keyt]=temp[key]

    #You don't want to flatten completely. You want to add the rows in one
    return grouped_GF_dic

def distancex(x,start_array,end_array):
    d=np.linalg.norm(x-start_array)/np.linalg.norm(end_array-start_array)
    return d

def grids_from_dic(stitched_dic,start,end):
    x_values=[]
    y_values=[]
    deg_values=[]
    start_karray=np.array([start[0],start[1]])
    end_karray=np.array([end[0],end[1]])
    def distance_one(x):
        return distancex(x,start_array=start_karray,end_array=end_karray)
    x_values=list(stitched_dic.keys())
    
    sorted_xvalues=sorted(x_values,key=distance_one)
    len_xvalues=[distance_one(x) for x in sorted_xvalues]
    for x in sorted_xvalues:
        y_values.append(stitched_dic[x][0])
        deg_values.append(stitched_dic[x][1])
    
    # print(len(sorted_xvalues))
    # print(len(y_values[499]))
    # print(len(deg_values[499]))
    


        
    return len_xvalues,y_values,deg_values


def grids_from_dicrot(stitched_dic,start,end):
    start_karray=np.array([start[0],start[1],start[2]])
    end_karray=np.array([end[0],end[1],end[2]])
    zgrid=[]
    x_values=[key[0] for key in list(stitched_dic.keys())]
    y_values=[key[1] for key in list(stitched_dic.keys())]#If stitched dic goes over all high symmetry paths I think you might run into problems with this because you cross the same y multiple times?
    # x=(0,0,0)
    # y=-10
    # test_list=[list(stitched_dic.keys())[0]]
    # print((x,3) in test_list)
    
    unique_xvalues=list(set(x_values))
    def distance_one(x):
        return distancex(x,start_array=start_karray,end_array=end_karray)
    sorted_unique_xvalues=sorted(unique_xvalues,key=distance_one)
    def distance_y(y):
        return (y-min(y_values))/(max(y_values)-min(y_values))
    
    sorted_yvalues=sorted(y_values,key=distance_y)
    print(len(stitched_dic.keys()))
    grid_arguments=[]
    for x in sorted_unique_xvalues:
        xrow=[]
        for key in stitched_dic.keys():
            if x==key[0]:
                xrow.append((x,key[1]))
        grid_arguments.append(xrow)
    
    for xrow in grid_arguments:
        zrow=[]
        for key in xrow:
            zrow.append(stitched_dic[key])
        zgrid.append(zrow)
    

            



    print([len(zrow) for zrow in grid_arguments])
    exit()

        
    return sorted_unique_xvalues,sorted_yvalues,zgrid
# test_group=group_parts_dic(directory=directory3,start_string="Gamma",end_string="Z",kpoints='20',parts='1',omegas="1000",dstr=dstr3)
# # print(list(test_group.keys())[9999+1])
# tx,ty,tz=grids_from_dic(stitched_dic=test_group,limit=1,start_karray=(0,0,0))
# print(len(tx))
# print(len(ty))
# print(len(tz[0]))
# print(tz[0][0])
# exit()

linescolored=16
linesplotted=16
linesinlabel=16
def chained_path_plot_link(path_list,start_path_count,axs,single,kpoints,directory,dstr,deg_list,mu_shift,params,variable,theta):
    deg_dict={}

    for i in range(0,len(path_list)-1):
        list_temp=[path_list[i],path_list[i+1]]
        start_string=var_file_reversed_Moire[tuple(list_temp[0])]

        end_string=var_file_reversed_Moire[tuple(list_temp[1])]
        grouped_dic=group_parts_dic(directory=directory,dstr=dstr,start_string=start_string,end_string=end_string,kpoints=kpoints,params=params,variable=variable,theta=theta)
        
        x_values,y_values,deg_values=grids_from_dic(stitched_dic=grouped_dic,start=tuple(list_temp[0]),end=tuple(list_temp[1]))
        

        eig_dict_loop={}
        for k in range(0,len(y_values[0])):
            #eig_dict_loop[k]=([y[k]-mu_shift-sorted(y)[0] for y in y_values],[j[k] for j in deg_values])
            eig_dict_loop[k]=([y[k]-mu_shift for y in y_values],[j[k] for j in deg_values])#[y[k]-y[0] for y in y_values]
            deg_list=list(set(list(set(deg_list))+list(set([j[k] for j in deg_values]))))
        # for j in eig_dict_loop.keys():
        for j in eig_dict_loop.keys():#The end-point sets how many lines you plot. [4,5]
            # annot = axs[i].annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
            #             bbox=dict(boxstyle="round", fc="w"),
            #             arrowprops=dict(arrowstyle="->"))
            # annot.set_visible(False)        
            #colorsmap=[colors[x] for x in eig_dict_loop[j][1]]
            colorsmap=[]
            for x in range(0,len(eig_dict_loop[j][1])):
                if j<linescolored:#This is to set how many lines you colors #UHK!#1.4 for not subtracted GS
                        colorsmap.append(colors[eig_dict_loop[j][1][x]])
                else:
                    # # kval=x_values[x]
                    # if (x_values[x]<0.7 and x_values[x]>0.3 and eig_dict_loop[j][0][x]<3.5):#1.5 for not subtracted GS #or (x_values[x]<0.8 and x_values[x]>0.6 and eig_dict_loop[j][0][x]<1.5):
                    #     colorsmap.append(colors[eig_dict_loop[j][1][x]])
                    # else:
                    colorsmap.append('dimgray')
                
            # print(eig_dict_loop[j][0][0])
            # print(x_values[0])
            if single:
                axs.scatter(x_values,eig_dict_loop[j][0],s=0.3,c=colorsmap)#facecolor=degcmap(eig_dict_loop[j][1])
                #axs.set_title(f"{var_dict_reversed[tuple(list_temp[0])]} to {var_dict_reversed[tuple(list_temp[1])]}",fontsize=tickfontsize)
                axs.set_xticks([-0.1])
                axs.set_xticklabels([var_dict_reversed_Moire[tuple(list_temp[0])]],fontsize=tickfontsize)#+ '\n'+ var_greek_dict[tuple(list_temp[0])]
                axs[i+start_path_count].tick_params(axis='x', length=0)
                axs[i+start_path_count].tick_params(axis='y', labelsize=ylabelsize)
                # axs.set_ylim([0, 2])
                if i==len(path_list)-2:
                    if np.equal(path_list[i+1],A).all():
                        axs[i+start_path_count].set_xticks([-0.1,1.1])
                        axs[i+start_path_count].set_xticklabels([var_dict_reversed_Moire[tuple(list_temp[0])],var_dict_reversed_Moire[tuple(list_temp[1])]],fontsize=tickfontsize)#+ '\n'+ var_greek_dict[tuple(list_temp[0])], + '\n'+ var_greek_dict[tuple(list_temp[1])]
                    else:
                        axs[i+start_path_count].set_xticks([-0.1,1.2])
                        axs[i+start_path_count].set_xticklabels([var_dict_reversed_Moire[tuple(list_temp[0])]+'|',var_dict_reversed_Moire[tuple(list_temp[1])]+'|'],fontsize=tickfontsize)#+ '\n'+ var_greek_dict[tuple(list_temp[0])], + '\n'+ var_greek_dict[tuple(list_temp[1])]
            else:
                axs[i+start_path_count].scatter(x_values,eig_dict_loop[j][0],s=0.3,c=colorsmap)#facecolor=degcmap(eig_dict_loop[j][1])
                #axs[i+start_path_count].set_title(f"{var_dict_reversed[tuple(list_temp[0])]} to {var_dict_reversed[tuple(list_temp[1])]}",fontsize=tickfontsize)
                axs[i+start_path_count].set_xticks([0])
                # g=var_greek_dict[tuple(list_temp[0])]
                #endstr=r'${{}}\tiny{}$'.format(g)
                # endstr=r'\small{{{}}}$'.format(g)
                axs[i+start_path_count].set_xticklabels([var_dict_reversed_Moire[tuple(list_temp[0])]+'\n'+ var_greek_dict_Moire[tuple(list_temp[0])]],fontsize=15)#+ 
                axs[i+start_path_count].tick_params(axis='x', length=0)
                axs[i+start_path_count].tick_params(axis='y', labelsize=ylabelsize)
                #axs[i+start_path_count].set_ylim([-11-6, -6-6])
                if i==len(path_list)-2:
                    if np.equal(path_list[i+1],A).all():
                        # print(np.equal(path_list[i+1],A).all())
                        axs[i+start_path_count].set_xticks([0,1])
                        axs[i+start_path_count].set_xticklabels([var_dict_reversed_Moire[tuple(list_temp[0])]+'\n'+ var_greek_dict_Moire[tuple(list_temp[0])],var_dict_reversed_Moire[tuple(list_temp[1])]+ '\n'+ var_greek_dict_Moire[tuple(list_temp[1])]],fontsize=15)
                        axs[i+start_path_count].tick_params(axis='x', length=0)
                    else:
                        if len(path_list)==2:
                            axs[i+start_path_count].set_xticks([0,1])
                            axs[i+start_path_count-1].set_xticks([1])#get rid of the last tick
                            axs[i+start_path_count+1].set_xticks([])
                            axs[i+start_path_count].set_xticklabels([r'$Z|X$',r'$R|M$'],fontsize=tickfontsize)
                            #axs[i+start_path_count].set_xticklabels([r'$|$'+var_dict_reversed[tuple(list_temp[0])],var_dict_reversed[tuple(list_temp[1])]+r'$|$'],fontsize=tickfontsize)#+ '\n'+ var_greek_dict[tuple(list_temp[0])], + '\n'+ var_greek_dict[tuple(list_temp[1])]
                            axs[i+start_path_count].tick_params(axis='x', length=0)
                        # else:
                        #     axs[i+start_path_count].set_xticks([-0.1,0.78])
                        #     axs[i+start_path_count].set_xticklabels([var_dict_reversed[tuple(list_temp[0])],var_dict_reversed[tuple(list_temp[1])]],fontsize=tickfontsize)#+ '\n'+ var_greek_dict[tuple(list_temp[0])], + '\n'+ var_greek_dict[tuple(list_temp[1])]
                        #     axs[i+start_path_count].tick_params(axis='x', length=0)
                
                # axs[i].set_ylim([0, 2])
                # if i==len(path_list)-2:
                #     axs[i+start_path_count].set_xticks([0,0.7])
                #     axs[i+start_path_count].set_xticklabels([var_dict_reversed[tuple(list_temp[0])]+ '-'+ var_greek_dict[tuple(list_temp[0])],var_dict_reversed[tuple(list_temp[1])]+ '-'+ var_greek_dict[tuple(list_temp[1])]],fontsize=5)
    # axs[5].set_xticks([])
    # axs[3].set_xticks([0])
    return deg_list

def chained_path_plot(path_lists,kpoints,directory,dstr,mu_shift,params,variable,theta):
    flatten_paths = [item for sublist in path_lists for item in sublist]
    deg_list=[]
    if len(flatten_paths)>2:
        single=False
        fig, axs = plt.subplots(1,len(flatten_paths)-len(path_lists),sharey=True)
        start_path_count=0
        for path_list in path_lists:
            deg_list=chained_path_plot_link(path_list=path_list,axs=axs,start_path_count=start_path_count,kpoints=kpoints,single=single,directory=directory,dstr=dstr,deg_list=deg_list,mu_shift=mu_shift,params=params,variable=variable,theta=theta)
            start_path_count=start_path_count+len(path_list)-1
    else:
        single=True
        fig,axs=plt.subplots(1)
        start_path_count=0
        for path_list in path_lists:
            chained_path_plot_link(path_list=path_list,axs=axs,start_path_count=start_path_count,kpoints=kpoints,single=single,directory=directory,dstr=dstr,mu_shift=mu_shift,params=params,variable=variable,theta=theta)
            start_path_count=0
    
    
    axs[0].set_ylabel(r'$E_{1}(\mathbf{k})/t$',fontsize=labelfontsize)
    #axs[0].set_ylabel(r'$(E_{4}(\mathbf{{k}})-E_{4,GS}(\mathbf{{k}}))/t$',fontsize=labelfontsize)
    plt.subplots_adjust(wspace=.0)
    custom_lines=[]
    names=[]
    for i in colors.keys():
        if i<=linesinlabel:#This sets how many lines you have in the deg label
            custom_lines.append(Line2D([0], [0], color=colors[i], lw=4))
            names.append(i)
    leg = plt.figlegend(custom_lines,names, loc=(0.87, 0.11),title="Degeneracy",fontsize=20,title_fontsize=20)
    #leg = plt.figlegend(custom_lines,names, loc=(0.87, 0.11),title="Degeneracy",fontsize=20,title_fontsize=20)#0.464
    fig.subplots_adjust(right=0.85)
    angle=theta.split('theta')[-1]
    fig.suptitle(r'One particle $\theta=$'+f'{angle}'+r'$^{\circ}$'+',  '+r'$U_{HK}=$'+f'{params[-1]}', fontsize=16)
    
    plt.show()

if __name__ == "__main__":
    mu_shift1=1
    params='UHK10'
    theta='theta1.05'
    #chained_path_plot(path_lists=[[Gamma,X,M]],kpoints="1000",directory=d6,dstr=dstr6,mu_shift=mu_shift1)
    chained_path_plot(path_lists=[[KM,GammaM,KM]],kpoints="100",directory=directory,dstr=dstr,mu_shift=mu_shift1,params=params,variable='oneparticle',theta=theta)
    # chained_path_plot([[Gamma,X,M,Gamma,Z],[X,R],[M,A,R,Z,A]])