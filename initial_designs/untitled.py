
#results=[]
#n_runs=20

#for i in range(n_runs):
#    print("Simulating {} of {} cells...".format(i+1,n_runs))
#    gen=gillespie(s0,t_G1,params)
#    results.append(gen)

#results=np.array(results)
#print(results[0])
#print(results[0][...,0])

fig=plt.figure(figsize=(12,8))

ax1=fig.add_subplot(5,1,1)
ax2=fig.add_subplot(5,1,2)
ax3=fig.add_subplot(5,1,3)
ax4=fig.add_subplot(5,1,4)
ax5=fig.add_subplot(5,1,5)

t_G1 /= 60

# plot all runs at once (the .T transposes the matrix)
ax1.plot(t_G1, results[:][...,0].T,'k:')
ax2.plot(t_G1, results[:][...,1].T,'y:')
ax3.plot(t_G1, results[:][...,2].T,'r:')
ax4.plot(t_G1, results[:][...,3].T,'b:')
ax5.plot(t_G1, results[:][...,4].T,'g:')

auxin_mean = np.average(results[:][...,0], axis=0)
TIR1_mean = np.average(results[:][...,1], axis=0)
RFP_mean = np.average(results[:][...,2],axis=0)
BFP_mean = np.average(results[:][...,3], axis=0)
volume_mean = np.average(results[:][...,4], axis=0)

print(np.average(RFP_mean))


ax1.plot(t_G1, auxin_mean,'k',label="Auxin")
ax2.plot(t_G1, TIR1_mean,'k',label="TIR1_mRNA")
ax3.plot(t_G1, RFP_mean,'k',label="TIR1")
ax4.plot(t_G1, BFP_mean,'k',label="RFP_mRNA")
ax5.plot(t_G1, volume_mean,'k',label="RFP")

ax1.set_xlabel("Time, in minutes")
ax1.set_ylabel("Concentration in umol per litre")
ax1.legend()

ax2.set_xlabel("Time, in minutes")
ax2.set_ylabel("Concentration in umol per litre")
ax2.legend()

ax3.set_xlabel("Time, in minutes")
ax3.set_ylabel("Concentration in umol per litre")
ax3.legend()

ax4.set_xlabel("Time, in minutes")
ax4.set_ylabel("Concentration in umol per litre")
ax4.legend()

ax5.set_xlabel("Time, in minutes")
ax5.set_ylabel("Volume in litres")
ax5.legend()



#ax1=fig2.add_subplot(5,1,1)
#ax2=fig2.add_subplot(5,1,2)
#ax3=fig2.add_subplot(5,1,3)
#ax4=fig2.add_subplot(5,1,4)
#ax5=fig2.add_subplot(5,1,5)


#t_all = np.linspace(0,G1_length*n_generations+G1_length,G1_length*n_generations+G1_length)

#ax1.plot(t_all, ori_results[:][...,0].flatten().tolist(),'g',label="Auxin")
#ax2.plot(t_all, ori_results[:][...,2].flatten().tolist(),'b',label="TIR1")
#ax3.plot(t_all, ori_results[:][...,4].flatten().tolist(),'y',label="RFP")
#ax4.plot(t_all, ori_results[:][...,6].flatten().tolist(),'r',label="BFP")
#ax5.plot(t_all, ori_results[:][...,7].flatten().tolist(),'r',label="Volume")

#ax1.set_xlabel("Time, in minutes")
#ax1.set_ylabel("Concentration in umol per litre")
#ax1.legend()

#ax2.set_xlabel("Time, in minutes")
#ax2.set_ylabel("Concentration in umol per litre")
#ax2.legend()

#ax3.set_xlabel("Time, in minutes")
#ax3.set_ylabel("Concentration in umol per litre")
#ax3.legend()

#ax4.set_xlabel("Time, in minutes")
#ax4.set_ylabel("Concentration in umol per litre")
#ax4.legend()

#ax5.set_xlabel("Time, in minutes")
#ax5.set_ylabel("Volume in litres")
#ax5.legend()




plt.show()