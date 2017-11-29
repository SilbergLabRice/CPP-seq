import numpy
import matplotlib
import scipy.stats as linstats
matplotlib.use('Agg')
from matplotlib import pyplot as pp
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib import ticker
import matplotlib.patches as mpatches

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Segoe UI']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['ytick.labelsize'] = 10
matplotlib.rcParams['xtick.labelsize'] = 10
matplotlib.rcParams['axes.labelsize'] = 12

def find_nearest(array,value):
    idx = (numpy.abs(array-value)).argmin()
    return array[idx]



genes = ['TnAK']#, 'GsAK', 'BsAK', 'BgAK']

for sample in genes:
    print('Plotting: ' + sample + '\n')
    figWidth = 12 #12 for squarish, 22 for wide
    figHeight = 24
    
    fig = pp.figure(figsize=(figWidth, figHeight))
    rows = 11

    #Import Data
    #Residue Counts
    selected_file = 'Outputs/' + sample + '_selected_ProteinNTerminiCounts.csv'
    unselected_file = 'Outputs/' + sample + '_unselected_ProteinNTerminiCounts.csv'

    selected = numpy.genfromtxt(selected_file, delimiter = ',', names=True, dtype=None)
    unselected = numpy.genfromtxt(unselected_file, delimiter = ',', names=True, dtype=None)
    
    #Foldchange
    fold_file = 'Outputs/' + sample + '_FoldChange.csv'
    fold = numpy.genfromtxt(fold_file, delimiter = ',', names=True, dtype=None)
    
    #Stats 
    stats_chosen = 'FoldTest'
    if stats_chosen == 'FoldTest':
        stats_file = 'Outputs/' + sample + '_FoldTest.csv'
        stats = numpy.genfromtxt(stats_file, delimiter = ',', names=True, dtype=None)
    if stats_chosen == 'FisherTest':
        stats_file = 'Outputs/' + sample + '_FisherTest.csv'
        stats = numpy.genfromtxt(stats_file, delimiter = ',', names=True, dtype=None)
    
    
    #Split CSVs
    selected_index = []
    selected_F1_parallel = []
    selected_F2_parallel = []
    selected_F3_parallel = []

    selected_F1_antiparallel = []
    selected_F2_antiparallel = []
    selected_F3_antiparallel = []

    for i in selected:
        selected_index.append(i[0])
        selected_F1_parallel.append(i[1])
        selected_F2_parallel.append(i[2])
        selected_F3_parallel.append(i[3])

        selected_F1_antiparallel.append(i[4])
        selected_F2_antiparallel.append(i[5])
        selected_F3_antiparallel.append(i[6])


    unselected_index = []
    unselected_F1_parallel = []
    unselected_F2_parallel = []
    unselected_F3_parallel = []

    unselected_F1_antiparallel = []
    unselected_F2_antiparallel = []
    unselected_F3_antiparallel = []

    for i in unselected:
        unselected_index.append(i[0])
        unselected_F1_parallel.append(i[1])
        unselected_F2_parallel.append(i[2])
        unselected_F3_parallel.append(i[3])

        unselected_F1_antiparallel.append(i[4])
        unselected_F2_antiparallel.append(i[5])
        unselected_F3_antiparallel.append(i[6])

    max_us = max(max([unselected_F1_parallel, unselected_F2_parallel, unselected_F3_parallel]))
    max_s = max(max([selected_F1_parallel, selected_F2_parallel, selected_F3_parallel]))
    
    fold_index = []
    fold_F1_parallel = []
    fold_F2_parallel = []
    fold_F3_parallel = []

    fold_F1_antiparallel = []
    fold_F2_antiparallel = []
    fold_F3_antiparallel = []
 
    for i in fold:
        fold_index.append(i[0])
        fold_F1_parallel.append(i[1])
        fold_F2_parallel.append(i[2])
        fold_F3_parallel.append(i[3])

        fold_F1_antiparallel.append(i[4])
        fold_F2_antiparallel.append(i[5])
        fold_F3_antiparallel.append(i[6])

    stats_index = []
    stats_F1_parallel = []
    stats_F2_parallel = []
    stats_F3_parallel = []
    stats_F1_antiparallel = []
    stats_F2_antiparallel = []
    stats_F3_antiparallel = []

    for i in stats:
        stats_index.append(i[0])
        stats_F1_parallel.append(i[1])
        stats_F2_parallel.append(i[2])
        stats_F3_parallel.append(i[3])
        stats_F1_antiparallel.append(i[4])
        stats_F2_antiparallel.append(i[5])
        stats_F3_antiparallel.append(i[6])
        
    #Generate colormap for statistics 
    sig_range = numpy.linspace(-300, -2, 50)
    rgb = cm.get_cmap(pp.get_cmap('Greys'))(numpy.linspace(0.5, 1.0, 50))[numpy.newaxis,:,:3][0][::-1]
    cMap = ListedColormap(rgb)
            
    ######################
    ####Plotting Loops####
    ######################
    frames = [1,2,3] 
    count = 1
    for i in frames:
        
        #import data 
        sel = globals()['selected_F' + str(i) + '_parallel']
        unsel = globals()['unselected_F' + str(i) + '_parallel']
        sel_ap = globals()['selected_F' + str(i) + '_antiparallel']
        unsel_ap = globals()['unselected_F' + str(i) + '_antiparallel']     
        stats_frame_parallel =  globals()['stats_F' + str(i)+'_parallel']
        stats_frame_antiparallel =  globals()['stats_F' + str(i)+'_antiparallel']
        fold_frame_parallel = globals()['fold_F' + str(i) + '_parallel']
        fold_frame_antiparallel = globals()['fold_F' + str(i) + '_antiparallel']
           
        

        stat_colors = []
        line_widths = []
        edge_colors = []
        for stat_sample in stats_frame_parallel:
            if stat_sample > -2:
                stat_colors.append([1,0,0,0.4])
                line_widths.append(0)
                edge_colors.append([1,0,0,1])
            if stat_sample <= -2:
                stat_colors.append(rgb[int(numpy.where(sig_range == find_nearest(sig_range, stat_sample))[0])])
                line_widths.append(0)
                edge_colors.append(rgb[int(numpy.where(sig_range == find_nearest(sig_range, stat_sample))[0])])
                
        stat_colors_antiparallel = []
        line_widths_antiparallel = []
        edge_colors_antiparallel = []
        for stat_sample in stats_frame_antiparallel:
            if stat_sample > -2:
                stat_colors_antiparallel.append([1,0,0,0.4])
                line_widths_antiparallel.append(0)
                edge_colors_antiparallel.append([1,0,0,1])
            if stat_sample <= -2:
                stat_colors_antiparallel.append(rgb[int(numpy.where(sig_range == find_nearest(sig_range, stat_sample))[0])])
                line_widths_antiparallel.append(0)
                edge_colors_antiparallel.append(rgb[int(numpy.where(sig_range == find_nearest(sig_range, stat_sample))[0])])
       
        ###REMOVE Unobserved Sequences###
        remove_unobserved = True
        observation_threshold = 0
        if remove_unobserved == False:
            observed_US_P_ind = range(0, len(unselected_index))
            observed_S_P_ind = range(0, len(unselected_index))
            observed_US_AP_ind = range(0, len(unselected_index))
            observed_S_AP_ind = range(0, len(unselected_index))
            observed_US_S_P = range(0, len(unselected_index))
            observed_US_S_AP = range(0, len(unselected_index))
            observed_P_AP = range(0, len(unselected_index))
            
            unobserved_US_P_ind = range(0, len(unselected_index))
            unobserved_S_P_ind = range(0, len(unselected_index))
            unobserved_US_AP_ind = range(0, len(unselected_index))
            unobserved_S_AP_ind = range(0, len(unselected_index))
            unobserved_US_S_P = range(0, len(unselected_index))
            unobserved_US_S_AP = range(0, len(unselected_index))
            unobserved_P_AP  = range(0, len(unselected_index))

        if remove_unobserved == True:
            #classify and remove unobserved positions 
            observed_US_P_ind = [unselected_index[index] for index, value in enumerate(unsel) if value > observation_threshold]
            observed_S_P_ind = [unselected_index[index] for index, value in enumerate(sel) if value > observation_threshold]
            observed_US_AP_ind = [unselected_index[index] for index, value in enumerate(unsel_ap) if value > observation_threshold]
            observed_S_AP_ind = [unselected_index[index] for index, value in enumerate(sel_ap) if value > observation_threshold]

            #observed in at least one of the two libraries 
            observed_US_S_P = map(int,list(set(observed_US_P_ind) | set(observed_S_P_ind)))
            observed_US_S_AP = map(int,list(set(observed_US_AP_ind) | set(observed_S_AP_ind))) 
            #observed in either insertion directions in either libraries, these samples should work best for Fisher's test
            observed_P_AP = list(set(observed_US_S_P) & set(observed_US_S_AP))


            #shift 0.5 just to line up with the barplots
            unobserved_US_P_ind = [unselected_index[index] for index, value in enumerate(unsel) if value <= observation_threshold]
            unobserved_S_P_ind = [unselected_index[index] for index, value in enumerate(sel) if value <= observation_threshold]
            unobserved_US_AP_ind = [unselected_index[index] for index, value in enumerate(unsel_ap) if value <= observation_threshold]
            unobserved_S_AP_ind = [unselected_index[index] for index, value in enumerate(sel_ap) if value <= observation_threshold]

            #unobserved across both of the two libraries 
            unobserved_US_S_P = list(set(unobserved_US_P_ind) & set(unobserved_S_P_ind))
            unobserved_US_S_AP = list(set(unobserved_US_AP_ind) & set(unobserved_S_AP_ind))

            #unobserved in at least one insertion directions in either of the two libraries, causes stastical issues with Fisher's test
            unobserved_P_AP = list(set(unobserved_US_S_P) | set(unobserved_US_S_AP))

            print('sum unobserved in either direction: ' + str(len(unobserved_P_AP)))
            print('sum of observed/unobserved: ' + str(len(observed_US_P_ind+unobserved_US_P_ind)))
            print('sum of observed/unobserved: ' + str(len(observed_P_AP+unobserved_P_AP)))
            print('unobserved in either P or AP' + str(len(unobserved_P_AP)))


            fold_frame_parallel = [fold_frame_parallel[x-1] for x in observed_P_AP]
            fold_frame_antiparallel = [fold_frame_antiparallel[x-1] for x in observed_P_AP]
            stats_frame_parallel = [stats_frame_parallel[x-1] for x in observed_P_AP]
            stats_frame_antiparallel = [stats_frame_antiparallel[x-1] for x in observed_P_AP]
        
        ######################
        ##RAW count plotting##
        ######################
        
        ax = pp.subplot(rows,len(frames),count)
        ax.bar(unselected_index, unsel, color = ['#00ff66'], width = 1, linewidth = 0, label = 'Parallel')
        ax.scatter(unobserved_US_P_ind, [0.01*(max(unsel)*1.05)]*len(unobserved_US_P_ind), s=2, facecolors='none', edgecolors='k', lw=0.25, label = 'Unobserved, $\emptyset$ = ' + str(len(unobserved_US_P_ind)) + '/' + str(len(unselected_index)))
        
        ax.set_xlim([0, max(unselected_index)])
        ax.set_ylim([-0.1, max(unsel)*1.05])
        ax.tick_params(direction = 'out', top = 'off', right = 'off')
        ax.legend(loc = 'upper center', fontsize = 8)
        ax.set_title('Frame: ' + str(count))
        ax.set_ylabel("Unselected Counts")
        ax.set_xlabel('Residue Position')
        
        ax2 = pp.subplot(rows,len(frames),count+len(frames))
        ax2.bar(selected_index, sel, color = ['#00ff66'], width = 1, linewidth = 0, label = 'Parallel') #rgb [0.40,1,0.50,1]
        ax2.scatter(unobserved_S_P_ind, [0.01*(max(sel)*1.05)]*len(unobserved_S_P_ind),s=2, facecolors='none', edgecolors='k', lw=0.25, label = 'Unobserved, $\emptyset$ = ' + str(len(unobserved_S_P_ind)) + '/' + str(len(selected_index)))
        
        ax2.set_xlim([0, max(selected_index)])
        ax2.set_ylim([-0.1, max(sel)*1.05])
        ax2.tick_params(direction = 'out', top = 'off', right = 'off')
        ax2.legend(loc = 'upper center', fontsize = 8)
        ax2.set_ylabel("Selected Counts")
        ax2.set_xlabel('Residue Position')

        ax = pp.subplot(rows,len(frames),count+len(frames)*2)
        ax.bar(unselected_index, unsel_ap, color = ['#ff0000'], width = 1, linewidth = 0, label = 'Antiparallel')
        ax.scatter(unobserved_US_AP_ind, [0.01*(max(unsel_ap)*1.05)]*len(unobserved_US_AP_ind),s=2, facecolors='none', edgecolors='k', lw=0.25, label = 'Unobserved, $\emptyset$ = ' + str(len(unobserved_US_AP_ind)) + '/' + str(len(unselected_index)))
        ax.set_xlim([0, max(unselected_index)])
        ax.set_ylim([-0.1, max(unsel_ap)*1.05])
        ax.tick_params(direction = 'out', top = 'off', right = 'off')
        ax.legend(loc = 'upper center', fontsize = 8)
        ax.set_ylabel("Unselected Counts")
        ax.set_xlabel('Residue Position')

        ax2 = pp.subplot(rows,len(frames),count+len(frames)*3)
        ax2.bar(selected_index, sel_ap, color = ['#ff0000'], width = 1, linewidth = 0, label = 'Antiparallel')
        ax2.scatter(unobserved_S_AP_ind, [0.01*(max(sel_ap)*1.05)]*len(unobserved_S_AP_ind), s=2, facecolors='none', edgecolors='k', lw=0.25, label = 'Unobserved, $\emptyset$ = ' + str(len(unobserved_S_AP_ind)) + '/' + str(len(selected_index)))
        ax2.set_xlim([0, max(selected_index)])
        ax2.set_ylim([-0.1, max(sel_ap)*1.05])
        ax2.tick_params(direction = 'out', top = 'off', right = 'off')
        ax2.legend(loc = 'upper center', fontsize = 8)
        ax2.set_ylabel("Selected Counts")
        ax2.set_xlabel('Residue Position')
        
        
        ################################
        #######Count Histograms#########
        ################################
        ax3 = pp.subplot(rows, len(frames),count+len(frames)*4)
        bins = numpy.logspace(0,5,50)
        ax3.hist([unsel_ap,unsel], bins = bins, histtype = 'barstacked', color = [(1,0,0,1),(0,1,0.4,1)], label = 'Unselected - Antiparallel')
        #ax3.hist(unsel, bins = bins, histtype = 'barstacked', color = (0,1,0.4,1), label = 'Unselected - Parallel')
        ax3.set_xscale('symlog')
        ax3.set_ylabel("Unique variants (#)")
        ax3.set_xlabel('Degenerate reads (#)')
        ax3.set_ylim([0, 30])
        
        ymin, ymax = ax3.get_ylim()
        xmin, xmax = ax3.get_xlim()
        ax3.text(xmax*10, (ymax/2), ['$\emptyset$ = ' + str(len(unobserved_US_AP_ind)) + '/' + str(len(selected_index))][0], ha = 'center', va = 'center', fontsize = 8, color = (1,0,0,1))
        ax3.text(xmax*10, ((ymax/2))*1.15, ['$\emptyset$ = ' + str(len(unobserved_US_P_ind)) + '/' + str(len(selected_index))][0], ha = 'center', va = 'center', fontsize = 8, color = (0,1,0.4,1))
        ax3.legend(loc='center left', bbox_to_anchor=(1, 0.3) , fontsize=7)
        
        ax3 = pp.subplot(rows, len(frames),count+len(frames)*5)
        bins = numpy.logspace(0,5,50)
        ax3.hist([sel_ap, sel], bins = bins, histtype = 'barstacked', color = [(1,0,0,1),(0,1,0.4,1)], label = 'Selected - Antiparallel')
        #ax3.hist(sel, bins = bins, histtype = 'step', color = (0,1,0.4,1), label = 'Selected - Parallel')
        ax3.set_xscale('symlog')
        ax3.set_ylabel("Unique variants (#)")
        ax3.set_xlabel('Degenerate reads (#)')
        ax3.set_ylim([0, 30])
        ymin, ymax = ax3.get_ylim()
        xmin, xmax = ax3.get_xlim()
        ax3.text(xmax*10, (ymax/2), ['$\emptyset$ = ' + str(len(unobserved_S_AP_ind)) + '/' + str(len(selected_index))][0], ha = 'center', va = 'center', fontsize = 8, color = (1,0,0,1))
        ax3.text(xmax*10, ((ymax/2))*1.15, ['$\emptyset$ = ' + str(len(unobserved_S_P_ind)) + '/' + str(len(selected_index))][0], ha = 'center', va = 'center', fontsize = 8, color = (0,1,0.4,1))
        ax3.legend(loc='center left', bbox_to_anchor=(1, 0.3) , fontsize=7)

        ################################
        ##Library Comparison plotting###
        ################################
        ##Adjust counts to only look at observed_P_AP
        
        sel = [sel[x-1] for x in observed_P_AP]
        unsel = [unsel[x-1] for x in observed_P_AP]
        sel_ap = [sel_ap[x-1] for x in observed_P_AP]
        unsel_ap = [unsel_ap[x-1] for x in observed_P_AP]
        stat_colors = [stat_colors[x-1] for x in observed_P_AP]
        stat_colors_antiparallel = [stat_colors_antiparallel[x-1] for x in observed_P_AP]
        edge_colors = [edge_colors[x-1] for x in observed_P_AP]
        

        print('Minimum Selected P count: ' + str(min(sel)))
        print('Minimum Selected AP count: '+ str(min(unsel)))
        print('Minimum Unselected P count: '+ str(min(sel_ap)))
        print('Minimum Unselected AP count: '+ str(min(unsel_ap)))
        
        ##Antiparallel vs Parallel 
        ax3 = pp.subplot(rows, len(frames),count+len(frames)*6)
        unsel_fit =  linstats.linregress(unsel_ap, unsel)
        unsel_fit_fn = numpy.poly1d([unsel_fit[0], unsel_fit[1]])
        
        sel_fit = linstats.linregress(sel_ap, sel)
        sel_fit_fn = numpy.poly1d([sel_fit[0], sel_fit[1]]) 
        
        x = [1,2,3,4]
        y= [1,2,3,4]
        equality =  linstats.linregress(x, y)
        equal_fn = numpy.poly1d([equality[0], equality[1]]) 

        ax3.scatter(sel_ap, sel, color = stat_colors, edgecolors = edge_colors, label = 'Selected: ' + str(len(sel)))
        ax3.scatter(unsel_ap, unsel, facecolors='none', edgecolors='k', label = 'Unselected: ' + str(len(unsel)))
        
        #Plot fit lines
        #equality
        xmax = int(max([max_s, max_us])*2)
        ax3.plot(range(-1,xmax), equal_fn(range(-1,xmax)), '-', color = '#000000', ms = 0)
        #unselected
        ax3.plot(range(-1,xmax), unsel_fit_fn(range(-1,xmax)), '--', color = '#6600ff', ms = 0)
        #selected
        ax3.plot(range(-1,xmax), sel_fit_fn(range(-1,xmax)), '--', color = '#ff6600', ms = 0)

        ax3.text(max([max_s, max_us])*200, 1250, ['y = ' + str(round(unsel_fit[0],3)) + '*x + ' + str(round(unsel_fit[1],3)) + '\n R^2: ' + str(round(unsel_fit[2]**2,3))][0], ha = 'center', va = 'center', fontsize = 8, color = '#6600ff')
        ax3.text(max([max_s, max_us])*200, 10000, ['y = ' + str(round(sel_fit[0],3)) + '*x + ' + str(round(sel_fit[1],3)) + '\n R^2: ' + str(round(sel_fit[2]**2,3))][0], ha = 'center', va = 'center', fontsize = 8, color = '#ff6600')
        
        #ax3.scatter(unsel, sel, color = [1./250.,166./250.,17./250.], label = 'Functional')
        #ax3.scatter(unsel_ap, sel_ap, color = [0.69, 0.13, 0.13], label = 'Non-Functional')
        ax3.set_xscale('symlog')
        ax3.set_yscale('symlog')
        ax3.set_xlim([-0.5, xmax])
        ax3.set_ylim([-0.5, xmax])
        ax3.tick_params(direction = 'out', top = 'off', right = 'off')
        ax3.legend(loc='center left', bbox_to_anchor=(1.125, 0.3) , fontsize=7)
        ax3.set_ylabel("Parallel Counts")
        ax3.set_xlabel('Antiparallel Counts')
        
        
        ##Selected v Unselected
        ax3 = pp.subplot(rows, len(frames),count+len(frames)*7)
        ap_fit =  linstats.linregress(unsel_ap, sel_ap)
        ap_fit_fn = numpy.poly1d([ap_fit[0], ap_fit[1]]) 
        
        #sel_fit = numpy.polyfit(unsel,sel, 1)
        p_fit =  linstats.linregress(unsel, sel)
        p_fit_fn = numpy.poly1d([p_fit[0], p_fit[1]]) 

        ax3.scatter(unsel, sel, color = stat_colors, edgecolors = edge_colors, label = 'Parallel: ' + str(len(sel)))
        ax3.scatter(unsel_ap,sel_ap, facecolors='none', edgecolors='k', label = 'Antiparallel: ' + str(len(sel_ap)))

        #Plot fit lines
        xmax = int(max([max_s, max_us])*2)
        #equality
        ax3.plot(range(-1,xmax), equal_fn(range(-1,xmax)), '-', color = '#000000', ms = 0)
        #parallel
        ax3.plot(range(-1,xmax), p_fit_fn(range(-1,xmax)), '--', color = '#00ff66', ms = 0)
        #antiparallel
        ax3.plot(range(-1,xmax), ap_fit_fn(range(-1,xmax)), '--', color = '#ff0000', ms = 0)
        
        ax3.set_xscale('symlog')
        ax3.set_yscale('symlog')
        ax3.set_xlim([-0.5, xmax])
        ax3.set_ylim([-0.5, xmax])
        ax3.tick_params(direction = 'out', top = 'off', right = 'off')
        ax3.legend(loc='center left', bbox_to_anchor=(1.125, 0.3) , fontsize=7)
        ax3.set_ylabel("Selected Counts")
        ax3.set_xlabel('Unselected Counts')
        
        ax3.text(max([max_s, max_us])*200, 1250, ['y = ' + str(round(ap_fit[0],3)) + '*x + ' + str(round(ap_fit[1],3)) + '\n R^2: ' + str(round(ap_fit[2]**2,3))][0], ha = 'center', va = 'center', fontsize = 8, color = '#ff0000')
        ax3.text(max([max_s, max_us])*200, 10000, ['y = ' + str(round(p_fit[0],3)) + '*x + ' + str(round(p_fit[1],3)) + '\n R^2: ' + str(round(p_fit[2]**2,3))][0], ha = 'center', va = 'center', fontsize = 8, color = '#00ff66')

        ######################
        #Fold-change plotting#
        ######################
        #replace nans with -9 for plotting fold change
        fold_frame_parallel = [-9.5 if numpy.isnan(x) else x for x in fold_frame_parallel]
        fold_frame_antiparallel = [-9.5 if numpy.isnan(x) else x for x in fold_frame_antiparallel]    
        
        ax3 = pp.subplot(rows, len(frames),count+len(frames)*8)
        #ax3.bar(selected_index, fold_frame_antiparallel, color = 'r', width = 1, linewidth = 0, label = 'Antiparallel')
        
        ax3.bar(observed_P_AP, fold_frame_parallel, color = stat_colors, width = 1, linewidth = 0, label = 'Parallel')
        ax3.scatter(unobserved_US_S_P, [0]*len(unobserved_US_S_P), s=2, facecolors='none', edgecolors='k', lw=0.25, label = 'Unobserved, $\emptyset$ = ' + str(len(unobserved_US_S_P)) + '/' + str(len(selected_index)))
        #ax3.scatter(unobserved_US_S_AP, [0]*len(unobserved_US_S_AP), s=1,  color = [1, 0, 0, 0.4], label = 'Unobserved, n = ' + str(len(unobserved_US_S_AP)) + '/' + str(len(selected_index)))
        
        ax3.axhline(y=-7, xmin=0, xmax=max(selected_index), linewidth=8, color = [0, 0, 0, 0.4], label = 'No Longer Observed')
        
        #axis adjustments
        ax3.axhline(y=0, xmin=0, xmax=max(selected_index), linewidth=0.25, color = [0, 0, 0, 1])
        ax3.set_xlim([0, max(selected_index)+1])
        ax3.set_ylim([-7, 7])
        ax3.tick_params(direction = 'out', top = 'off', right = 'off')
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(2)) #tickspacing
        ax3.legend(loc='center left', bbox_to_anchor=(1.25, 0.5) , fontsize=7)
        
       #legend,swap range
        sm = pp.cm.ScalarMappable(cmap=cMap, norm=pp.Normalize(vmin=int(min(sig_range)), vmax=int(max(sig_range))))
        # fake up the array of the scalar mappable. Urgh...
        sm._A = []
        #formatter = ticker.LogFormatter(10, labelOnlyBase=False)
        cb = pp.colorbar(sm)#, format = formatter)
        font_size = 7 # Adjust as appropriate.
        cb.ax.tick_params(labelsize=font_size)
        tick_locator = ticker.MaxNLocator(nbins=7)
        cb.locator = tick_locator
        cb.update_ticks()
        
        
        ax3.set_ylabel("Fold Change (log2)")
        ax3.set_xlabel('Residue Position')

        ax3 = pp.subplot(rows, len(frames),count+len(frames)*9)
        ax3.bar(observed_P_AP, fold_frame_antiparallel, color = stat_colors_antiparallel, width = 1, linewidth = 0, label = 'Antiparallel')
        ax3.scatter(unobserved_US_S_AP, [0]*len(unobserved_US_S_AP), s=2, facecolors='none', edgecolors='k', lw=0.25, label = 'Unobserved, $\emptyset$ = ' + str(len(unobserved_US_S_AP)) + '/' + str(len(selected_index)))
        
        ax3.axhline(y=-7, xmin=0, xmax=max(selected_index), linewidth=8, color = [0, 0, 0, 0.4], label = 'No Longer Observed')
        
        #axis adjustments
        ax3.axhline(y=0, xmin=0, xmax=max(selected_index), linewidth=0.25, color = [0, 0, 0, 1])
        ax3.set_xlim([0, max(selected_index)+1])
        ax3.set_ylim([-7, 7])
        ax3.tick_params(direction = 'out', top = 'off', right = 'off')
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(2)) #tickspacing
        ax3.legend(loc='center left', bbox_to_anchor=(1.25, 0.5) , fontsize=7)
        
        ax3.set_ylabel("Fold Change (log2)")
        ax3.set_xlabel('Residue Position')
        
        
        #legend,swap range
        sm = pp.cm.ScalarMappable(cmap=cMap, norm=pp.Normalize(vmin=int(min(sig_range)), vmax=int(max(sig_range))))
        # fake up the array of the scalar mappable. Urgh...
        sm._A = []
        #formatter = ticker.LogFormatter(10, labelOnlyBase=False)
        cb = pp.colorbar(sm)#, format = formatter)
        font_size = 7 # Adjust as appropriate.
        cb.ax.tick_params(labelsize=font_size)
        tick_locator = ticker.MaxNLocator(nbins=7)
        cb.locator = tick_locator
        cb.update_ticks()
        

        
        ###################################
        #Fold-change v. Abundance plotting#
        ###################################
        #unsel = [unsel[x-1] for x in observed_US_S_P]
        #unsel_ap = [unsel_ap[x-1] for x in observed_US_S_P]
        
        #Average Unselected and Selected counts for x-axis
        #use Unselected counts for x-axis
        #remove no longer observed 'nan' to calculate mean/stdev
        background = [x for x in fold_frame_antiparallel if x > -9] 
        background_mean = numpy.average(background)
        print('mean: ' + str(background_mean))
        background_stdev = numpy.std(background)
         
        greaterThanBackground = [fold_frame_parallel[ind] for ind, x in enumerate(fold_frame_parallel) if x >= background_mean + 2*background_stdev]
        greaterThanBackgroundx = [unsel[ind] for ind, x in enumerate(fold_frame_parallel) if x >= background_mean + 2*background_stdev]
        
        ax3 = pp.subplot(rows, len(frames),count+len(frames)*10)
        
        ax3.scatter(unsel, fold_frame_parallel, color = stat_colors, edgecolors = edge_colors, label = 'Parallel')
        ax3.scatter(unsel_ap, fold_frame_antiparallel, facecolors='none', edgecolors='k', label = 'Antiparallel')

        print('# of parallel above 2sigma: ' + str(len(greaterThanBackground)) + '\n')
        
        xmin, xmax = ax3.get_xlim()
        ax3.axhline(y=background_mean, xmin=xmin, xmax=xmax, linewidth=1, color = 'k', label = 'mean')
        ax3.axhline(y=-9.5, xmin=xmin, xmax=xmax, linewidth=8, color = [0, 0, 0, 0.4], label = 'No Longer Observed')

        ax3.legend(loc='center left', bbox_to_anchor=(1.05, 0.5) , fontsize=7)
        ax3.set_ylim([-10, 6])
        ax3.tick_params(direction = 'out', top = 'off', right = 'off')
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(2)) #tickspacing
        ax3.set_xlim([-0.5, max([max_us]) * 1.5])
        ax3.set_xscale('symlog')
        ax3.set_ylabel("Fold Change (log2)")
        ax3.set_xlabel('Unselected Counts')
        
        
        count = count + 1
        
    pp.tight_layout()
    fig.subplots_adjust(wspace=1)
    figname = 'Plots/' + sample + '_Correlation Plot.svg'
    pp.savefig(figname, bbox_inches="tight")



    