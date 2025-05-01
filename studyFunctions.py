#Estimate cell volume based on Beauvoit et al
def estimateVcell(T):
    import math
    Vcell = ((0.02894333/(0.00760074+(math.e**(-0.18324543*T))))-0.03277816)*(10**-12) # units = m3
    return Vcell

#Estimate pericarp volume based on Jose's data
def estimateVpericarp(T,hollow=False):
    #assuming sphere
    Vpericarp = (0.0000889527 +((0.000003593996 - 0.0000889527)/(1 +((T/28.56023)**9.692893)))) # units = m3
    #
    if hollow:
        #assuming hollowsphere - 20%
        Vpericarp = (0.0000007121165 +((0.00000002819613 - 7.121165e-7)/(1 +((T/28.55635)**9.567273)))) # units = m3
    return Vpericarp

#Estimate number of cells
def estimateNcell(T,hollow=False):
    Vcell = estimateVcell(T)				#unit = m3
    Vpericarp = estimateVpericarp(T,hollow)	#unit = m3
    return Vpericarp/Vcell

#Estimate volume of cell wall in a cell
def estimateVcellwall(T):
    Vcell = estimateVcell(T)			#unit = m3
    a = Vcell**(1.0/3)
    #cell wall thickness = 100 nm according to figure 19-68, https://www.ncbi.nlm.nih.gov/books/NBK26928/
    b=100*(10**-9)				#unit = m
    y = (8*(b**3))+(6*(a**2)*b)-(12*a*(b**2))
    return y


#Estimate volume of cell wall
def estimateTotalVcellwall(T,Ncell=0,hollow=False):
    Vcellwall = estimateVcellwall(T)		#unit = m3
    if Ncell == 0:
        Ncell = estimateNcell(T,hollow)
    return Vcellwall*Ncell			#unit = m3/fruit

#Esimate cellulose content
def celluloseLevels(T,Ncells=0,hollow=False):
    p_cellulose = 1.54*(10**6)              # unit = g/m3
    x = estimateTotalVcellwall(T,Ncells,hollow)    # unit = m3/fruit
    y = x*p_cellulose                       # unit = g/fruit
    z = (y/180.16) * 1000                   # unit = mmol/fruit
    return z


#Estimate cellulose demand flux
def celluloseDemandFlux(t,Ncells=0,unit_time=1):
    z = celluloseLevels(t,Ncells)-celluloseLevels(t-unit_time,Ncells)	#unit = mmol/fruit/unit_time
    return z

#Estimate cell membrane of a cell
def estimateVcellmembrane(T):
    Vcell = estimateVcell(T)
    a=Vcell**(1.0/3)
    #cell wall thickness = 100 nm according to figure 19-68, https://www.ncbi.nlm.nih.gov/books/NBK26928/
    b=100*(10**-9)
    #cell membrane thickness = 10 nm according to https://hypertextbook.com/facts/2001/JenniferShloming.shtml (see Hine, Robert. "Membrane." The Facts on File Dictionary of Biology. 3rd ed. New York: Checkmark, 1999: 198.)
    #c = 10*(10**-9)				#unit = m
    #cell membrane thickness = 4 nm BioNumbers
    c = 4*(10**-9)				#unit = m
    y = (8*(c**3))+(6*(a**2)*c)+(24*(b**2)*c)-(12*(c**2)*(a-(2*b)))-(24*a*b*c)	#unit = m3
    return y

#Estimate total volume of cell membrane
def estimateTotalVcellmembrane(T,Ncell=0,hollow=False):
    Vcellmembrane = estimateVcellmembrane(T)		#unit = m3
    if Ncell == 0:
        Ncell = estimateNcell(T,hollow)
    return Vcellmembrane*Ncell			#unit = m3/fruit

#Esimate phospholipid content
def phospholipidLevels(T,Ncells=0,scaling_factor=1):
    p_phospholipid = 0.975*(10**6)              	# unit = g/m3
    x = estimateTotalVcellmembrane(T,Ncells)    	# unit = m3/fruit
    y = x*p_phospholipid                       		# unit = g/fruit
    z = (y/689.978637415773) * 1000 * scaling_factor	# unit = mmol/fruit
    return z

#Estimate phospholipid demand flux
def phospholipidDemandFlux(t,Ncells=0,unit_time=1,scaling_factor=1):
    z = phospholipidLevels(t,Ncells,scaling_factor)-phospholipidLevels(t-unit_time,Ncells,scaling_factor)	#unit = mmol/fruit/unit_time
    return z


#Estimate total cytosolic volume
def estimateVcyt(T,Ncell=0,hollow=False):
    import math
    j = T*24*60
    Vvac = 0.853*(1-math.e**(float(-2292-j)/10633))
    Vcyt = (0.933-Vvac)/1.13
    if Ncell == 0:
        Ncell = estimateNcell(T,hollow)
    return Vcyt*estimateVcell(T)*Ncell

#Estimate total protein content
def estimateProteinContent(T,Ncell=0,ProtConc=21458.1747597,hollow=False):
    if Ncell == 0:
        Ncell = estimateNcell(T,hollow)
    Vcyt = estimateVcyt(T,Ncell,hollow)		# unit = m3/fruit
    #ProtConc#                                  # unit = mmol/m3
    return ProtConc*Vcyt                        # unit = mmol/fruit

#Estimate Protein demand
def estimateProteinDemandFlux(T,Ncell=0,ProtConc=21458.1747597,unit_time=1,hollow=False):
    if Ncell == 0:
        Ncell = estimateNcell(T,hollow)
    z = estimateProteinContent(T,Ncell,ProtConc=21458.1747597) - estimateProteinContent(T-unit_time,Ncell,ProtConc=21458.1747597)	# unit = mmol/fruit/unit_time
    return z

#Estimate fruit C content
def estimateCcontent(T,hollow=False):
    y = (0.03617391*(estimateVpericarp(T,hollow)/estimateVpericarp(59,hollow))*100)+0.18173913
    return y

#Estimate phloem uptake rate
def estimatePhloemUptakeConstraint(Ccont):
    a=1.21863077
    b=6.70786808
    c=0.43646842
    y = a+b/(Ccont+c)
    return y			# mgC/fruit/hr
