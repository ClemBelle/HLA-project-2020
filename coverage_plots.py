## IMPORT
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import sys
import xlsxwriter

## GLOBAL VARIABLES
# Class I - Exons position
HLAA=[[29942532,29942626],[29942757,29943026],[29943268,29943543],[29944122,29944397],[29944500,29944616],[29945059,29945090],[29945234,29945281],[29945451,29945870]]
HLAB=[[31357086,31357179],[31356687,31356957],[31356167,31356442],[31355317,31355592],[31355107,31355223],[31354633,31354665],[31354483,31354526],[31353875,31354296]]
HLAC=[[31271999,31272086],[31271599,31271868],[31271073,31271347],[31270210,31270485],[31269966,31270085],[31269493,31269525],[31269338,31269385],[31268749,31269173]]
# Class II - Exons position
HLADQA1=[[32637406,32637540],[32641310,32641558],[32641972,32642253],[32642610,32642784],[32642952,32643671]]
HLADQB1=[[32666499,32666657],[32664798,32665067],[32661967,32662247],[32661347,32661457],[32660859,32660882],[32659467,32660249]]
HLADRB1=[[32589643,32589836],[32584109,32584378],[32581557,32581838],[32580746,32580856],[32580247,32580270],[32578769,32579104]]

HLA_TOTAL=[HLAA,HLAB,HLAC,HLADQA1,HLADQB1,HLADRB1]

# HLA names
Names=['A','A_filt','B','B_filt','C','C_filt','DQA1','DQA1_filt','DQB1','DQB1_filt','DRB1','DRB1_filt']
Segments=['exon 1','intron 1','exon 2', 'intron 2', 'exon 3', 'intron 3', 'exon 4', 'intron 4', 'exon 5', 'intron 5', 'exon 6', 'intron 6', 'exon 7', 'intron 7', 'exon 8']

# Exons color
Colors=[['pink','orange','lime','palegreen','lightseagreen','skyblue','violet','chocolate','lightgrey'],['deeppink','orangered','limegreen','mediumseagreen','teal','dodgerblue','darkmagenta','maroon','darkgrey']]
# First sub-list is for the "low" quality reads, so the ones obtained by making the difference between total and filtered reads
# Second sub-list is for the "high" quality reads, the ones directly obtained with samtools depth by using a MAPQ filter

## DATA & CALCULS
def get_content(sample_ID):
    """Retrieve the content of the 12 coverage files"""
    TOTAL=[]
    file_name=sample_ID+'-cov_depth_HLA-'
    for i in Names:
        file_name=sample_ID+'-cov_depth_HLA-'+i+".coverage"
        file=open(file_name,"r")
        CONTENT=[]
        for i in file.readlines():
            CONTENT=CONTENT+[i.split()]
        file.close()
        CONTENT=CONTENT[1:]                 #Remove the header
        for i in CONTENT:
            i[1],i[2]=int(i[1]),int(i[2])   #Convert the position and nb of reads in integers
        TOTAL=TOTAL+[CONTENT]
    return TOTAL


def mean_inter(L,inter):
    """Groups the bases together and calculates the means per group"""
    X=[]
    for k in range(len(L)):
        Xsub=[]
        Lind=[]
        for i in range(0,len(L[k]),inter):
            Lind=Lind+[i]
        for i in Lind[:-1]:
            mean=0
            for j in range(i,i+inter):
                mean=mean+L[k][j][2]
            mean=mean/inter
            Xsub=Xsub+[round(mean,2)]
        mean=0
        c=0
        for i in range(Lind[-1],len(L)):
            mean=mean+L[k][i][2]
            c=c+1
        if c!=0:
            mean=mean/c
            Xsub=Xsub+[round(mean,2)]
        X=X+[Xsub]
    
    return X
    
def correction_total(X):
    """Does the difference between all reads and filtered reads. 
    Replace the total read count by that difference, so the total will be the addition of those 2"""
    for i in range(1,len(X),2):
        for j in range (len(X[i])):
            X[i-1][j]=round(X[i-1][j]-X[i][j],2)
    return X


def base_color(L,HLA_TOTAL,Lcolor):
    """Attributes a color per base. 
    One color is attributed to a specific exon or introns"""
    Lbase=[]
    
    for i in range(len(L)):             # Go through the different HLA
        if i%2==0:
            ind=int(i/2)
        else:
            ind=int((i-1)/2)
        
        Lbasesub=[]
        for j in range(len(L[i])):      # Go through the content of an HLA 
            found=0
            for k in range(len(HLA_TOTAL[ind])):   # Go through the exons coordinates of that HLA gene
                if L[i][j][1]>=HLA_TOTAL[ind][k][0] and L[i][j][1]<=HLA_TOTAL[ind][k][1]:
                    Lbasesub=Lbasesub+[Lcolor[i%2][k]]     #Lcolor is in the same order as the exon, i%2=0 for all the non filtered and =1 for the filtered, colors are in the right order
                    found=1
            if found==0:            # Means the coordinates don't belong to an exon
                Lbasesub=Lbasesub+[Lcolor[i%2][-1]]    #Exon color is the last of Lcolor
        Lbase=Lbase+[Lbasesub]
    return Lbase
    
def exon_mean(L, Lbase):
    """"Calculates the mean per segments : mean per exon and mean per intron
    Considers a segment finished at each color change"""
    Lmean=[]
    for i in range(len(Lbase)):
        Lmeansub=[]
        mean=L[i][0][2]
        count_ex=1
        count=0
        for j in range(1,len(Lbase[i])):
            if Lbase[i][j]==Lbase[i][j-1]:
                mean=mean+L[i][j][2]
                count_ex=count_ex+1
            else:
                count=count+count_ex
                Lmeansub=Lmeansub+[round(mean/count_ex,2)]
                mean=L[i][j][2]
                count_ex=1
        # For the last exon : 
        last_ex=len(Lbase[i])-count
        mean=0
        for k in range(-last_ex,0):
            mean=mean+L[i][k][2]
        Lmeansub=Lmeansub+[round(mean/last_ex,2)]   
        Lmean=Lmean+[Lmeansub]
    for j in [2,3,4,5,8,9,10,11]: #Reverse HLA genes
        Lmean[j].reverse()
    return Lmean
    
def color_inter(Lbase,inter,Lcolor):
    """Attributes a color to each columns of the bar plot. 
    If one column contains intron + exon, it will be colored with the exon's color"""
    
    Lcol=[]
    for i in range(len(Lbase)): #each HLA list
        Lcolsub=[]
        Lind=[]
        for j in range(0,len(Lbase[i]),inter):
            Lind=Lind+[j]
        
        for k in Lind[:-1]:
            color=Lcolor[i%2][-1]
            for h in range(k,k+inter):
                if Lbase[i][h]!=Lcolor[i%2][-1]:
                    color=Lbase[i][h]
            Lcolsub=Lcolsub+[color]
        
        color=Lcolor[i%2][-1]
        for l in range(Lind[-1],len(Lbase[i])):
            if Lbase[i][l]!=Lcolor[i%2][-1]:
                color=Lbase[i][l]
        Lcolsub=Lcolsub+[color]
            
        Lcol=Lcol+[Lcolsub]
    return Lcol
        
def coord_mean_exon(HLA,inter):
    """Retrieves coordinates to draw the mean line over the exons' columns on the bar plot"""
    Coord=[]
    if HLA[0][0]<HLA[-1][0]:
        ref=HLA[0][0]
    else:
        ref=HLA[-1][0]
    for i in range(len(HLA)):
        Coord=Coord+[[round(((HLA[i][0]-ref)/inter),2),round(((HLA[i][1]-ref)/inter),2)]]
    
    for j in Coord:
        if round(j[0],0)>j[0]:
            j[0]=round(j[0],0)-1.45
        elif round(j[0],0)<=j[0]:
            j[0]=round(j[0],0)-0.45
        if round(j[1],0)>j[1]:
            j[1]=round(j[1],0)-0.55
        elif round(j[1],0)<=j[1]:
            j[1]=round(j[1],0)+0.45
        
    return Coord
        
def create_array(Lm,Lc, HLA_TOTAL,ind,inter):
    """Transforms the coordinates of the mean lines into an array to be plotted"""
    
    c=0
    coord=[]
    for i in range(0,len(Lm[ind*2]),2): #Only exons coord
        coord1=[Lc[c][0],Lm[ind*2][i]]
        coord2=[Lc[c][1],Lm[ind*2][i]]
        coord=coord+[[coord1]+[coord2]]
        c=c+1
    return coord

    
def coord_annotation(Coord, Lm,ind,inter):
    Lann=[]
    for i in range(0,len(Lm[ind*2]),2):
        if i%2==0:
            x=Coord[int(i/2)][0]
        else:
            x=Coord[int((i-1)/2)][0]
        y=Lm[ind*2][i]
        Lann=Lann+[[x,y]]
    return Lann
    
## PLOT


def bar_plot(Xtot,Xfilt,Lcolor,Lcolorf,Lmean,HLA_TOTAL,ind,sample_ID,HLA,qual,inter):
    """Creates the plot"""
    Title=sample_ID+": HLA-"+HLA+" coverage"
    # Total height 
    bars = np.add(Xfilt, Xtot).tolist()
    # Create the bottom bars (high quality ones)
    plt.bar(range(len(Xtot)), Xfilt, color=Lcolorf, width=0.9)
    # Create the top bars ("low quality" ones)
    plt.bar(range(len(Xtot)), Xtot, bottom=Xfilt, color=Lcolor, width=0.9)
    # Add a grid in the background
    plt.grid(axis='y', alpha=0.5)
    
    # Label axis and title
    plt.ylabel('Nb of reads',weight='bold')
    plt.xlabel('*'+str(inter)+' pb\n\nHigh quality = MAPQ > '+qual)
    plt.title(Title,weight='bold')
    
    #Always there
    intron_bar = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'darkgrey')
    intron_bar_filt = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'lightgrey')
    exon1_bar = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'deeppink')
    exon1_bar_filt = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'pink')
    exon2_bar = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'orangered')
    exon2_bar_filt = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'orange')
    exon3_bar = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'limegreen')
    exon3_bar_filt = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'lime')
    exon4_bar = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'mediumseagreen')
    exon4_bar_filt = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'palegreen')
    exon5_bar = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'teal')
    exon5_bar_filt = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'lightseagreen')
    if ind==3:
        plt.legend([intron_bar, intron_bar_filt, exon1_bar, exon1_bar_filt, exon2_bar, exon2_bar_filt, exon3_bar, exon3_bar_filt, exon4_bar, exon4_bar_filt, exon5_bar, exon5_bar_filt], ['intron high quality', 'intron', 'exon 1 high quality','exon 1 ', 'exon 2high quality','exon 2', 'exon 3 high quality', 'exon 3', 'exon 4 high quality','exon 4', 'exon 5 high quality','exon 5'], markerscale = 100, frameon = True, fontsize = 5)
        
    elif ind in [0,1,2,4,5]:
        exon6_bar = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'dodgerblue')
        exon6_bar_filt = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'skyblue')
        if ind in [4,5]:
            plt.legend([intron_bar, intron_bar_filt, exon1_bar, exon1_bar_filt, exon2_bar, exon2_bar_filt, exon3_bar, exon3_bar_filt, exon4_bar, exon4_bar_filt, exon5_bar, exon5_bar_filt, exon6_bar, exon6_bar_filt], ['intron high quality', 'intron', 'exon 1 high quality','exon 1 ', 'exon 2high quality','exon 2', 'exon 3 high quality', 'exon 3', 'exon 4 high quality','exon 4', 'exon 5 high quality','exon 5', 'exon 6 high quality','exon 6'], markerscale = 100, frameon = True, fontsize = 5)
            
        elif ind in [0,1,2]:
            exon7_bar = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'darkmagenta')
            exon7_bar_filt = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'violet')
            exon8_bar = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'crimson')
            exon8_bar_filt = matplotlib.patches.Rectangle((0, 0), 0, 0, color = 'lightcoral')
            plt.legend([intron_bar, intron_bar_filt, exon1_bar, exon1_bar_filt, exon2_bar, exon2_bar_filt, exon3_bar, exon3_bar_filt, exon4_bar, exon4_bar_filt, exon5_bar, exon5_bar_filt, exon6_bar, exon6_bar_filt, exon7_bar, exon7_bar_filt, exon8_bar, exon8_bar_filt ], ['intron high quality', 'intron', 'exon 1 high quality','exon 1 ', 'exon 2high quality','exon 2', 'exon 3 high quality', 'exon 3', 'exon 4 high quality','exon 4', 'exon 5 high quality','exon 5', 'exon 6 high quality','exon 6', 'exon 7 high quality','exon 7', 'exon 8 high quality','exon 8'], markerscale = 100, frameon = True, fontsize = 5)
        
    Lc=coord_mean_exon(HLA_TOTAL[ind],inter)
    lines = np.array(create_array(Lmean,Lc, HLA_TOTAL,ind,inter))
    x, y = lines.T
    plt.plot(x, y, linewidth=1, linestyle='--',color='black')
    Lann=coord_annotation(Lc, Lmean,ind,inter)
    if ind in [2,3,4,5]:
        Lann.reverse()
    for k in Lann:
        plt.annotate(str(k[1]),(k[0],(k[1]+1)),fontsize=7)
    
    
        
    

def plot_assemble(sample_ID, inter, qual, X, Lcol, Lmean, HLA_TOTAL): 
    """Assembles all 6 graphs together and generates a PDF as output"""
    
    Name_PDF=sample_ID+'_HLA_coverage_plots.pdf'
    
    f=plt.figure(figsize=(20,17), dpi=100)
    
    plt.subplot(3,2,1)
    bar_plot(X[0],X[1],Lcol[0],Lcol[1],Lmean, HLA_TOTAL,0,sample_ID,"A",qual,inter) 
    
    plt.subplot(3,2,2)
    bar_plot(X[6],X[7],Lcol[6],Lcol[7],Lmean, HLA_TOTAL,3,sample_ID,"DQA1",qual,inter)
    
    plt.subplot(3,2,3)
    bar_plot(X[2],X[3],Lcol[2],Lcol[3],Lmean, HLA_TOTAL,1,sample_ID,"B",qual,inter)
    
    plt.subplot(3,2,4)
    bar_plot(X[8],X[9],Lcol[8],Lcol[9],Lmean, HLA_TOTAL,4,sample_ID,"DQB1",qual,inter)
    
    plt.subplot(3,2,5)
    bar_plot(X[4],X[5],Lcol[4],Lcol[5],Lmean, HLA_TOTAL,2,sample_ID,"C",qual,inter)
    
    plt.subplot(3,2,6)
    bar_plot(X[10],X[11],Lcol[10],Lcol[11],Lmean, HLA_TOTAL,5,sample_ID,"DRB1",qual,inter)
    
    plt.subplots_adjust(hspace=0.5)
    plt.show()
    
    f.savefig(Name_PDF,bbox_inches='tight')
    
## CSV 
def output_csv_mean(Lmean,Names,Segments,sample_ID):
    """Creates a .csv file to store the mean coverage of each segments of the HLA genes"""
    name=sample_ID+"_coverage_means_HLA_genes.xlsx"
    workbook = xlsxwriter.Workbook(name)
    worksheet = workbook.add_worksheet()
    bold = workbook.add_format({'bold': True})
    
    # Header
    row = 0
    col = 1
    for HLA in Names:
        worksheet.write(row, col,"HLA-"+HLA,bold)
        col+=1
        
    #1st column
    row=1
    col=0
    for segment in Segments:
        worksheet.write(row,col,segment,bold)
        row+=1
    
    #Data
    row=1
    col=1
    for i in Lmean: 
        for mean in i: 
            worksheet.write(row,col,mean)
            row+=1
        col+=1
        row=1
    workbook.close()
    

## MAIN 

#Questions for the user
sample_ID=sys.argv[1]   #str(input("Sample ID : "))
inter=30                #int(input("How many bases do you want to regroup ? : ")) #30 is a good one
qual='50'               #str(input('What cut-off did you use to filter your data ? : '))


# Retrieve the content from the files
TOTAL=get_content(sample_ID)

# Group some bases together + Calculates the mean reads per group + Nb of reads filtered, nb of reads in total
X=mean_inter(TOTAL,inter)
X=correction_total(X)

# Create the colors lists
Lbase=base_color(TOTAL,HLA_TOTAL,Colors)

# Calculates the mean per segments (exons/introns)
Lmean=exon_mean(TOTAL,Lbase)

#Attributes a color per column
Lcol_int=color_inter(Lbase,inter,Colors)

#Plots the data
plot_assemble(sample_ID, inter, qual, X, Lcol_int, Lmean, HLA_TOTAL)

#Stores the means in a xlsx file
output_csv_mean(Lmean,Names,Segments,sample_ID)

