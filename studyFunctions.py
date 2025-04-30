#This function processes raw biomass data and returns upper and lower bounds to
#be used in constraint scanning as done by Colombie et al 2015
#
#
def generateBoundsFromBiomass(datafile="/home/sanu/ColobieDataRaw.csv",met="sucrose",Nsampling = 1000,DPA=(11.4,14.8,18.2,21.6,25,28.4,31.8,25.2,38.6,41.3),show_plots=True,start=4.0,stop=57.0,degree=3,Ssampling=0.75):
  import pandas as pd
  import math
  import numpy as np

  conc=list()
  conc_max = list()
  conc_min = list()
  rate_max=list()
  rate_min=list()

  #import data
  df = pd.read_csv(datafile,sep="\t")
  num2rem = int(0.025*Nsampling)

  x_values = list()
  y_values = list()


  for i in range(0,len(df)):
    if not df[met][i]== 0:
      x_values.append(df["DPA"][i])
      y_values.append(df[met][i])


  log_y_values = list()
  for i in y_values:
    log_y_values.append(math.log(i))



  x2 = np.arange(start,stop+0.1,0.1)
  y2 = np.poly1d(np.polyfit(x_values,log_y_values,degree))

  ys = dict()

  for i in range(0,Nsampling):
    ind = np.random.choice(range(0,len(x_values)),size=int(len(x_values)*Ssampling),replace=False)
    #print ind
    x=list()
    y=list()
    for j in ind:
      x.append(x_values[j])
      y.append(log_y_values[j])
    ys[i]=np.poly1d(np.polyfit(x,y,degree))


  maxys = list()
  max95s = list()
  minys = list()
  min95s = list()
  temp=dict()
  for i in ys.keys():
    temp[i]=ys[i](x2)

  yi = dict()
  for j in range(0,len(x2)):
    templist=list()
    for i in temp.keys():
      templist.append(temp[i][j])
    temp2list=sorted(templist)[num2rem:Nsampling-num2rem]
    yi[j]=templist
    maxys.append(max(templist))
    max95s.append(max(temp2list))
    minys.append(min(templist))
    min95s.append(min(temp2list))




  y3 = list()
  for i in y2(x2):
    y3.append(math.exp(i))

  maxy1 = list()
  for i in maxys:
    maxy1.append(math.exp(i))


  miny1 = list()
  for i in minys:
    miny1.append(math.exp(i))

  maxy95_1 = list()
  for i in max95s:
    maxy95_1.append(math.exp(i))


  miny95_1 = list()
  for i in min95s:
    miny95_1.append(math.exp(i))

  #print y3
  #print "----"
  #print x2

  for i in DPA:
    #print i
    conc.append(y3[int((i-min(x_values))*10)])
    conc_max.append(maxy95_1[int((i-min(x_values))*10)])
    conc_min.append(miny95_1[int((i-min(x_values))*10)])

  if show_plots:
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 20}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=3 # makes axes line thicker
    plt.figure(figsize=(5,5))


    ax = plt.subplot()
    ax.plot(x_values,y_values,".",label="raw data")
    ax.plot(x2,y3,"g-",label="fitted curve")
    ax.plot(x2,maxy1,"r--",label="fitted curve with max. y")
    ax.plot(x2,miny1,"b--",label="fitted curve with min. y")
    plt.xlabel("Time(days)")
    plt.ylabel("Concentration (micromol/fruit)")
    plt.legend(bbox_to_anchor=(2, 1),fontsize=15)
    plt.show()

  ##############################derivatives#####################

  y2_deriv = y2.deriv()


  ys_deriv = dict()
  for i in ys.keys():
    ys_deriv[i] = ys[i].deriv()

  maxys_deriv = list()
  maxy95_deriv = list()
  minys_deriv = list()
  miny95_deriv = list()
  temp=dict()
  for i in ys.keys():
    templist=list()
    for j in range(0,len(x2)):
      templist.append(y3[j]*ys_deriv[i](x2[j]))
    temp[i]=templist


  yi_deriv = dict()
  for j in range(0,len(x2)):
    templist=list()
    for i in temp.keys():
      templist.append(temp[i][j])
    temp2list=sorted(templist)[num2rem:Nsampling-num2rem]
    yi_deriv[j]=templist
    maxys_deriv.append(max(templist))
    maxy95_deriv.append(max(temp2list))
    minys_deriv.append(min(templist))
    miny95_deriv.append(min(temp2list))




  y4 = list()
  for i in range(0,len(x2)):
    y4.append(y3[i]*y2_deriv(x2)[i])

  rate_max = list()
  rate_min = list()
  if not len(DPA)==0:
    j=0
    for i in DPA:
      rate_max.append(maxy95_deriv[int((i-min(x_values))*10)])
      rate_min.append(miny95_deriv[int((i-min(x_values))*10)])
      j=j+1


  if show_plots:
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 20}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=3 # makes axes line thicker
    plt.figure(figsize=(5,5))


    ax = plt.subplot()
    ax.plot(x2,y4,"g-",label="derivitative of fitted curve")
    ax.plot(x2,maxys_deriv,"r-",label="derivative of fitted curve with max.y")
    ax.plot(x2,minys_deriv,"b-",label="derivative of fitted curve with min.y")
    ax.plot(x2,maxy95_deriv,"r--",label="derivative of fitted curve with max. y 95%")
    ax.plot(x2,miny95_deriv,"b--",label="derivative of fitted curve with min. y 95%")
    ax.set_title(met+" flux")
    plt.xlabel("Time(days)")
    plt.ylabel("Flux (micromol/fruit/day)")
    plt.legend(bbox_to_anchor=(2.55, 1),fontsize=15)
    plt.show()

  return (conc,conc_max,conc_min,rate_max,rate_min)


#Remove copy tag
def removeCopyTag(ID):
  ID=ID.replace("10","").replace("2","").replace("3","").replace("4","").replace("5","").replace("6","").replace("7","").replace("8","").replace("9","").replace("1","")
  return ID


#Function to plot biomass accumulation rates as barplots with upper and lower bo-
#-unds as the upper and lower ends of the y-axis
#
#args: 1) a list of reaction objects
#output: a plot
def plotMetabAccumulationRates(rxnlist):
  cols = 6
  import matplotlib.pyplot as plt
  plt.rcParams.update({'font.size': 20}) #sets a global fontsize
  plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
  plt.rcParams['xtick.major.width'] = 1
  plt.rcParams['ytick.major.size'] = 5
  plt.rcParams['ytick.major.width'] = 1
  plt.rcParams['axes.linewidth']=3 # makes axes line thicker
  plt.figure(figsize=(20,10))

  rows = int(len(rxnlist)/cols)
  if len(rxnlist)%cols != 0:
    rows=rows+1

  i=0
  for rxn in rxnlist:
    i=i+1
    ax = plt.subplot(rows,cols,i)
    ax.bar(1,rxn.x,0.5)
    ax.set_ylim([rxn.lower_bound,rxn.upper_bound])
    ax.set_xlim([0.75,1.75])
    ax.set_title(rxn.name.replace("Biomass_",""))
    ax.set_xticks([])
    ax.axhline(y=0,color="grey")

  plt.tight_layout()
  plt.show()


#this function converts a model which fractionally charged metabolites to a classical model
def convertToClassicalModel(cobra_model2,comp="*",updateCharges=""):
    if not updateCharges == "":
        #override charges from file
        fin = open(updateCharges,"r")
        ChargeDict=dict()
        for line in fin:
            met=line.replace("\n","").split("\t")[0]
            met = met.replace("-","_")
            charge = line.replace("\n","").split("\t")[1]
            ChargeDict[met]=charge
        fin.close()

        temp=cobra_model2.copy()
        for met in temp.metabolites:
            if(met.compartment=="b"):
                cobra_model2.metabolites.get_by_id(met.id).remove_from_model()

        for met in cobra_model2.metabolites:
            tempMet=met.id
            if(["1","2","3","4","5","6","7","8","9"].__contains__(met.id[len(met.id)-1])):
                tempMet = met.id[0:len(met.id)-1]
            elif(met.id[len(met.id)-2:]=="10"):
                tempMet = met.id[0:len(met.id)-2]
            if(ChargeDict.keys().__contains__(tempMet)):
                met.charge = ChargeDict.get(tempMet)
            if met.charge is None:
                met.charge=0

    import re
    #List metabolites that were added
    fractionMets=set()
    for rxn in cobra_model2.reactions:
        for met in rxn.metabolites.keys():
            a=re.search("^a{1,3}",met.id)
            anion=""
            if a:
                anion=a.group(0)
            b=re.search("^b{1,3}",met.id)
            basic=""
            if b:
                basic=b.group(0)
            if (abs(rxn.metabolites.get(met)) % 1 > 0 and (not anion == "" or not basic== "")):
                if not comp == "*":
                    if met.compartment == comp:
                        fractionMets.add(met.id)
                else:
                    fractionMets.add(met.id)
    cobra_model = cobra_model2.copy()

    for met in fractionMets:
        for rxn in cobra_model.metabolites.get_by_id(met).reactions:
            #print rxn.id
            a=re.search("^a{1,3}",met)
            anion=""
            if a:
                anion=a.group(0)
            b=re.search("^b{1,3}",met)
            basic=""
            if b:
                basic=b.group(0)
            prefix = anion
            if prefix == "":
                prefix = basic
            #print(met)
            #if(rxn.id=="Protein_Processing_c1"):
                #print rxn.reaction

            coeff1 = rxn.metabolites.get(cobra_model.metabolites.get_by_id(met[len(prefix):]))
            rxn.add_metabolites({cobra_model.metabolites.get_by_id(met[len(prefix):]):-1*coeff1})
            coeff2 = rxn.metabolites.get(cobra_model.metabolites.get_by_id(met))
            rxn.add_metabolites({cobra_model.metabolites.get_by_id(met):-1*coeff2})
            Charge=(coeff1*float(cobra_model.metabolites.get_by_id(met[len(prefix):]).charge))+(coeff2*float(cobra_model.metabolites.get_by_id(met).charge))-((coeff1+coeff2)*float(cobra_model.metabolites.get_by_id(met[len(prefix):]).charge))
            if rxn.id.__contains__("_Transfer"):
                compSet = set()
                for m in rxn.metabolites:
                    compSet.add(m.compartment)
                this = cobra_model.metabolites.get_by_id(met).compartment
                #print this
                rxn.add_metabolites({cobra_model.metabolites.get_by_id(met[len(prefix):]):coeff1+coeff2,cobra_model.metabolites.get_by_id("PROTON_"+cobra_model.metabolites.get_by_id(met).compartment):round(Charge,4)})
                other = compSet.difference(this)
                for c in other:
                    #print c
                    if cobra_model.metabolites.get_by_id("PROTON_"+c).reactions.__contains__(rxn):
                        prot_bal = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+c)) + rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+this))
                        if prot_bal > 0:
                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+c))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+c):-1*tempCoeff})

                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+this))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+this):-1*tempCoeff})

                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+c):prot_bal})
                        elif prot_bal < 0:
                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+c))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+c):-1*tempCoeff})

                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+this))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+this):-1*tempCoeff})

                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+c):prot_bal})
                        else:
                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+c))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+c):-1*tempCoeff})

                            tempCoeff = rxn.metabolites.get(cobra_model.metabolites.get_by_id("PROTON_"+this))
                            rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_"+this):-1*tempCoeff})
            else:
                rxn.add_metabolites({cobra_model.metabolites.get_by_id(met[len(prefix):]):coeff1+coeff2,cobra_model.metabolites.get_by_id("PROTON_"+cobra_model.metabolites.get_by_id(met).compartment):round(Charge,4)})
            #if(rxn.id=="Protein_Processing_c1"):
            #    print Charge
            #    print rxn.reaction
    uncharged = cobra_model.copy()
    return uncharged

def convertToFractionalChargesFruit(uncharged,infile="MetaboliteChargedStates.xlsx",compH={"v1":5.5}):
    model = uncharged.copy()

    #from xlrd import open_workbook
    import openpyxl
    wb = openpyxl.load_workbook(infile)  
    sheet = wb.active
    from cobra.core import Metabolite

    #wb = open_workbook(infile)

    FractionDict = dict()
    FractionCharge = dict()
    for comp in compH.keys():
        pH = compH.get(comp)
        tempDict1 = dict()
        for sheet in wb.sheets():
            met = str(sheet.name)
            metID = met+"_"+comp
            tempDict2 = dict()
            for i in range(1,sheet.ncols):
                charge = sheet.cell(0,i).value
                row = int((pH*10)-9)
                pH = sheet.cell(row,0).value
                n = int(model.metabolites.get_by_id(metID).charge) - charge
                prefix = ""
                if n < 0:
                    prefix = "a"*int(abs(n))
                elif n > 0:
                    prefix = "b"*int(abs(n))
                met_spcl = prefix+met
                FractionCharge[met_spcl]=int(charge)
                fra = sheet.cell(row,i).value
                if fra > 0:
                    if tempDict2.keys().__contains__(met_spcl):
                        tempDict2[met_spcl] = tempDict2[met_spcl]+fra
                    else:
                        tempDict2[met_spcl] = fra
            tempDict1[met]=tempDict2
        FractionDict[pH] = tempDict1

    rxnUpdated = set()

    for comp in compH.keys():
        pH = compH.get(comp)
        for met in FractionDict[pH].keys():
            for rxn in model.metabolites.get_by_id(met+"_"+comp).reactions:
                if rxn.id.__contains__("Transfer"):
                    rxnUpdated.add(rxn)
                coeff = rxn.metabolites.get(model.metabolites.get_by_id(met+"_"+comp))
                rxn.add_metabolites({model.metabolites.get_by_id(met+"_"+comp):-1*coeff})
                for met_spcl in FractionDict[pH][met].keys():
                    coeff_new = coeff*FractionDict[pH][met][met_spcl]*0.01
                    #print(met_spcl)
                    #print(coeff_new)
                    f = 0
                    for m in model.metabolites.query(met_spcl+"_"+comp):
                        #print m
                        if m.id==met_spcl+"_"+comp:
                            f=1
                            break
                    if f == 0:
                        m=Metabolite(met_spcl+"_"+comp)
                        m.name = model.metabolites.get_by_id(met+"_"+comp).name
                        m.formula = model.metabolites.get_by_id(met+"_"+comp).formula
                        if "a" in met_spcl:
                            if m.formula.__contains__("H"):
                                parts = m.formula.split("H")
                                m.formula = parts[0]+"H"+str(int(parts[1][0])+(len(met_spcl)))+parts[1][1:]
                            else:
                                m.formula = m.formula+"H"+(len(met_spcl)*1)
                        else:
                            if m.formula.__contains__("H"):
                                parts = m.formula.split("H")
                                m.formula = parts[0]+"H"+str(int(parts[1][0])-(len(met_spcl)))+parts[1][1:]
                        m.notes = model.metabolites.get_by_id(met+"_"+comp).notes
                        m.charge = FractionCharge[met_spcl]
                        m.compartment = comp
                    rxn.add_metabolites({m:coeff_new})
                Charge = 0
                for tmet in rxn.metabolites.keys():
                    Charge = Charge + (int(tmet.charge) * rxn.metabolites.get(tmet))
                if Charge != 0:
                    #print(rxn.reaction)
                    #print Charge
                    rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+model.metabolites.get_by_id(met+"_"+comp).compartment):-1*Charge})
                    #print(rxn.reaction)

    for rxn in rxnUpdated:
        #print(rxn.reaction)
        compSet = set()
        n = 0
        for met in rxn.metabolites.keys():
            if met.id.__contains__("PROTON"):
                n=n+1
            compSet.add(met.compartment)
        #if n < 2:
        #    continue
        c=list(compSet)
        a=c[0][len(c[0])-1]
        b=c[1][len(c[1])-1]
        if a > b:
            thisC = c[1]
            nextC = c[0]
        else:
            thisC = c[0]
            nextC = c[1]
        #if rxn.id.__contains__("MAL_v_Transfer"):
        #    print("This ="+thisC)
        #    print("Next ="+nextC)
        #    print(rxn.reaction)
        if rxn.metabolites.keys().__contains__(model.metabolites.get_by_id("PROTON_"+thisC)):
            this_prot = rxn.metabolites.get(model.metabolites.get_by_id("PROTON_"+thisC))
            tempCoeff = rxn.metabolites.get(model.metabolites.get_by_id("PROTON_"+thisC))
            rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+thisC):-1*tempCoeff})
        else:
            this_prot = 0
        if rxn.metabolites.keys().__contains__(model.metabolites.get_by_id("PROTON_"+nextC)):
            next_prot = rxn.metabolites.get(model.metabolites.get_by_id("PROTON_"+nextC))
            tempCoeff = rxn.metabolites.get(model.metabolites.get_by_id("PROTON_"+nextC))
            rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+nextC):-1*tempCoeff})
        else:
            next_prot = 0
        net = this_prot+next_prot
        if net > 0:
            rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+nextC):net})
        elif net < 0:
            rxn.add_metabolites({model.metabolites.get_by_id("PROTON_"+thisC):net})
        #if rxn.id.__contains__("MAL_v_Transfer"):
        #    print(rxn.reaction)
    return model

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

#Estimate fruit FW
def estimateFruitFW(T):
    a=60.81973868579338
    b=-10657.387676309034
    c=-10598.32971442436
    d=32.132292282374394
    e=8.381671975168445
    y =a+((b-c)/(1+((T/d)**e)))
    return y

#Estimate Lycopene content
def estimateLycopene(T):
    a=8.02771901
    b=31961.5923
    c=31969.6384
    d=44.0003073
    e=18.0700518
    y = (a +((b - c)/(1 +((T/d)**e)))) # units = m3
    return y

#Estimate Phytol content
def estimatePhytol(DPA):
    if DPA<18:
        z = 0.004
    elif DPA>=18 and DPA<38:
        z = (DPA*0.06804829)-1.09814513
    elif DPA>=38 and DPA<42:
        z = (DPA*(-0.29816473))+12.56436756
    else:
        z = 0
    return z


#Estimate Lycopene demand
def estimateLycopeneDemand(DPA,unit_time):
    return estimateLycopene(DPA)-estimateLycopene(DPA-unit_time)

#Estimate Phytol demand
def estimatePhytolDemand(DPA,unit_time):
    return estimatePhytol(DPA) - estimatePhytol(DPA-unit_time)

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

#Function to convert any model with reversible reactions to a copy of the same m-
#-odel with only irreversible reactions. ID of reverse reactions are generated by
#suffixing "_reverse" to the ID of the orignal reaction.
#args: 1) a cobra model
#output: a cobra model with only irreversible reactions
def rev2irrev(cobra_model):
  exp_model=cobra_model.copy()

  for RXN in cobra_model.reactions:
    rxn=exp_model.reactions.get_by_id(RXN.id)
    if (rxn.lower_bound < 0):
      rxn_reverse = rxn.copy()
      rxn_reverse.id = "%s_reverse" %(rxn.id)
      rxn.lower_bound = 0
      rxn_reverse.upper_bound = 0
      exp_model.add_reaction(rxn_reverse)

  return exp_model



#Function to constraint sum of fluxes when performing FBA
#args: 1) a cobra model, 2) a python list of reactions to leave out from constrai-
#-nt, 3) the float value that sum of fluxes must be constrained to & 4) value obj-
#-ective function needs to be constraint to (provide "" to avoid constraining obj-
#ective function)
#output: a cobra model with sum of fluxes constrained to
def constrainSumOfFluxes(cobra_model, rxn2avoid,SFvalue,objvalue):
  temp=cobra_model.copy()
  SFMet = Metabolite("SFMet",name="Sum of fluxes pseudometabolite",compartment="c2")
  for rxn in cobra_model.reactions:
    if not rxn2avoid.__contains__(rxn.id):
        if rxn.id.__contains__("reverse"):
            temp.reactions.get_by_id(rxn.id).add_metabolites({SFMet:-1})
        else:
          temp.reactions.get_by_id(rxn.id).add_metabolites({SFMet:1})
    SFRxn = Reaction("SFRxn",name="Sum of fluxes pseudoreaction")
    SFRxn.add_metabolites({SFMet:-1})
    SFRxn.lower_bound=SFvalue
    SFRxn.upper_bound=SFvalue
    temp.add_reaction(SFRxn)
    if (not objvalue=="") and (len(temp.objective) == 1):
        for rxn in temp.objective.keys():
          rxn.lower_bound=objvalue
          rxn.upper_bound=objvalue
  return temp


#################################################################################
# This function is a modified version of cobrapy pfba function			#
#										#
#################################################################################

import logging
from warnings import warn
from itertools import chain

from optlang.symbolics import Zero

from cobra.util import solver as sutil
from cobra.core.solution import get_solution

def pfba_Weighted(model, weightings, fraction_of_optimum=1.0, objective=None, reactions=None):
    """Perform basic pFBA (parsimonious Enzyme Usage Flux Balance Analysis)
    to minimize total flux.
    pFBA [1] adds the minimization of all fluxes the the objective of the
    model. This approach is motivated by the idea that high fluxes have a
    higher enzyme turn-over and that since producing enzymes is costly,
    the cell will try to minimize overall flux while still maximizing the
    original objective function, e.g. the growth rate.
    Parameters
    ----------
    model : cobra.Model
        The model
    fraction_of_optimum : float, optional
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum.
    objective : dict or model.problem.Objective
        A desired objective to use during optimization in addition to the
        pFBA objective. Dictionaries (reaction as key, coefficient as value)
        can be used for linear objectives.
    reactions : iterable
        List of reactions or reaction identifiers. Implies `return_frame` to
        be true. Only return fluxes for the given reactions. Faster than
        fetching all fluxes if only a few are needed.
    Returns
    -------
    cobra.Solution
        The solution object to the optimized model with pFBA constraints added.
    References
    ----------
    .. [1] Lewis, N. E., Hixson, K. K., Conrad, T. M., Lerman, J. A.,
       Charusanti, P., Polpitiya, A. D., Palsson, B. O. (2010). Omic data
       from evolved E. coli are consistent with computed optimal growth from
       genome-scale models. Molecular Systems Biology, 6,
       390. doi:10.1038/msb.2010.47
    """
    reactions = model.reactions if reactions is None \
        else model.reactions.get_by_any(reactions)
    with model as m:
        add_pfba_Weighted(m, weightings, objective=objective,
                 fraction_of_optimum=fraction_of_optimum)
        m.slim_optimize(error_value=None)
        solution = get_solution(m, reactions=reactions)
    return solution


#################################################################################
# This function is a modified version of cobrapy add_pfba function			#
#										#
#################################################################################

def add_pfba_Weighted(model, weightings, objective=None, fraction_of_optimum=1.0):
    """Add pFBA objective
    Add objective to minimize the summed flux of all reactions to the
    current objective.
    See Also
    -------
    pfba
    Parameters
    ----------
    model : cobra.Model
        The model to add the objective to
    objective :
        An objective to set in combination with the pFBA objective.
    fraction_of_optimum : float
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum.
    """
    if objective is not None:
        model.objective = objective
    if model.solver.objective.name == '_pfba_objective':
        raise ValueError('The model already has a pFBA objective.')
    sutil.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
    reaction_variables = ((rxn.forward_variable, rxn.reverse_variable)
                          for rxn in model.reactions)
    variables = chain(*reaction_variables)
    model.objective = model.problem.Objective(
        Zero, direction='min', sloppy=True, name="_pfba_objective")
    #print([v for v in variables])
    tempDict = dict()
    for v in variables:
        w = str(v).split("=")[1].replace(" ","").replace("<","")
        found=False
        for rxn in weightings.keys():
            if w.__contains__(rxn):
                #print(v)
                #print(rxn)
                tempDict[v]=weightings[rxn]
                found=True
                break
        if not found:
            print("Weightings for reaction "+w+" not found, so assuming weighting = 1")
            tempDict[v] = 1
    model.objective.set_linear_coefficients(tempDict)



#Function to perform FVA analysis which maintains sum of fluxes at a minimal val-
#ue
#args: 1) a cobra model 2) Objective 3) reaction to avoid when constraining sum
#of fluxes 4) reaction list for FVA 5) solver used to perform FVA
#output: a cobra model with FVA as an attribute called fva
def FBA_FVA_run(cobra_model,obj,rxn2avoid = [],rxnlist=[],solver="",weightings={}):
  from cobra import flux_analysis
  if len(weightings)>0:
    weightings_submittied=True
  else:
    weightings_submittied=False

  if len(rxnlist)==0:
    rxnlist = cobra_model.reactions
  for rxn in cobra_model.reactions:
    if not rxn.id in weightings.keys():
      weightings[rxn.id]=1
      if weightings_submittied:
          print("Warning")
          print(rxn.id)
          return
  print("Runing pFBA")
  solution = pfba_Weighted(cobra_model,weightings)
  objvalue = solution.fluxes.get(obj.id)
  a = 0
  for i in cobra_model.reactions:
    a = a + abs(solution.fluxes.get(i.id)*weightings[i.id])

  sumOfFluxes = a

  cobra_model2 = cobra_model.copy()
  irr_model = rev2irrev(cobra_model2)
  print("Setting SOF model")
  sfmodel = constrainSumOfFluxes(irr_model,rxn2avoid,sumOfFluxes,objvalue,weightings)
  rxnlist2 = list()
  if rxnlist == cobra_model.reactions:
    rxnlist2 = sfmodel.reactions
  else:
    for rxn in rxnlist:
      if rxn.lower_bound<0 and rxn.upper_bound>0:
        rxnlist2.append(sfmodel.reactions.get_by_id(rxn.id+"_reverse"))
      rxnlist2.append(sfmodel.reactions.get_by_id(rxn.id))
  #print("Rxn list ="+str(rxnlist2))
  tempsol = sfmodel.optimize()
  print(tempsol.fluxes.get(obj.id))
  print("Running FVA")

  if solver != "":
    import optlang
    if optlang.available_solvers.keys().__contains__(solver) and optlang.available_solvers[solver]:
      sfmodel.solver=solver
    else:
      print("Requested solver "+solver+" not available, using current model solver...")
  fva = flux_analysis.flux_variability_analysis(sfmodel,reaction_list = rxnlist2)
  print("Processing results")

  fva2=dict()
  for mode in fva.keys():
    if mode == "maximum":
      tempdict = dict()
      FVArxnSet = set()
      for rxn in fva[mode].keys():
        if rxn.__contains__("_reverse"):
          rxn = rxn.replace("_reverse","")
        if FVArxnSet.__contains__(rxn):
          continue
        FVArxnSet.add(rxn)
        if not fva[mode].keys().__contains__(rxn+"_reverse"):
          maxi = fva[mode][rxn]
        else:
          maxi = fva[mode][rxn]+fva[mode][rxn+"_reverse"]
        tempdict[rxn]=maxi
    else:
      tempdict=dict()
      FVArxnSet = set()
      for rxn in fva[mode].keys():
        if rxn.__contains__("_reverse"):
          rxn = rxn.replace("_reverse","")
        if FVArxnSet.__contains__(rxn):
          continue
        FVArxnSet.add(rxn)
        if not fva[mode].keys().__contains__(rxn+"_reverse"):
          mini = fva[mode][rxn]
        else:
          mini = fva[mode][rxn]+fva[mode][rxn+"_reverse"]
        tempdict[rxn]=mini
    fva2[mode]=tempdict

  sfmodel.fva = fva
  cobra_model.fva = fva2
  cobra_model.solution = solution
  return cobra_model


#Function to convert any model with reversible reactions to a copy of the same m-
#-odel with only irreversible reactions. ID of reverse reactions are generated by
#suffixing "_reverse" to the ID of the orignal reaction.
#args: 1) a cobra model
#output: a cobra model with only irreversible reactions
def rev2irrev(cobra_model):
  exp_model=cobra_model.copy()

  for RXN in cobra_model.reactions:
    rxn=exp_model.reactions.get_by_id(RXN.id)
    if (rxn.lower_bound < 0):
      rxn_reverse = rxn.copy()
      rxn_reverse.id = "%s_reverse" %(rxn.id)
      rxn.lower_bound = 0
      rxn_reverse.upper_bound = 0
      exp_model.add_reaction(rxn_reverse)

  return exp_model


#Function to constraint sum of fluxes when performing FBA
#args: 1) a cobra model, 2) a python list of reactions to leave out from constrai-
#-nt, 3) the float value that sum of fluxes must be constrained to & 4) value obj-
#-ective function needs to be constraint to (provide "" to avoid constraining obj-
#ective function) 5) Flux weightings
#output: a cobra model with sum of fluxes constrained to
def constrainSumOfFluxes(cobra_model, rxn2avoid,SFvalue,objvalue,weightings):
  from cobra.core import Metabolite, Reaction

  temp=cobra_model.copy()
  SFMet = Metabolite("SFMet",name="Sum of fluxes pseudometabolite",compartment="c2")
  for rxn in cobra_model.reactions:
    if not rxn2avoid.__contains__(rxn.id):
        if rxn.id.__contains__("reverse"):
            temp.reactions.get_by_id(rxn.id).add_metabolites({SFMet:-1*weightings[rxn.id.replace("_reverse","")]})
        else:
            temp.reactions.get_by_id(rxn.id).add_metabolites({SFMet:1*weightings[rxn.id.replace("_reverse","")]})
    SFRxn = Reaction("SFRxn",name="Sum of fluxes pseudoreaction")
    SFRxn.add_metabolites({SFMet:-1})
    SFRxn.lower_bound=SFvalue
    SFRxn.upper_bound=SFvalue
    temp.add_reaction(SFRxn)
    if (not objvalue=="") and (len([rxn for rxn in temp.reactions if rxn.objective_coefficient==1]) == 1):
        for rxn in [rxn for rxn in temp.reactions if rxn.objective_coefficient==1]:
          rxn.lower_bound=float(objvalue)
          rxn.upper_bound=float(objvalue)
  return temp


#####################################################################
# This function generates ATP budgets for a given flux distribution #
# inputs: 1) an FBA model, 2) a dictionary object with reaction ids #
# as keys and reaction fluxes as values, 3) name of output file (op-#
# -tional), 4) Option to show plots, 5) If choosing to show plot, c-#
# -hoose wether to use percentage or absolute values in the plot. 6)#
# Provide a day or night indicator tag to specify day or night ATP  #
# summary 7) a destination file to save plot to 8) a dictionary to  #
# specify colour for fluxes in plot                                 #
#####################################################################
def generateATPbudget(model,solution,outfile="",show_plot=True,percentage=False,day_or_night_tag="1",save_plot_to="temp.png",colourDict={}):
  if outfile!="":
    fout = open(outfile,"w")
  ATPdict = dict()
  total = 0
  for p in ("c","p","m","x"):
    met=model.metabolites.get_by_id("ATP_"+p+day_or_night_tag)
    met1=model.metabolites.get_by_id("aATP_"+p+day_or_night_tag)
    for rxn in met.reactions:
      if rxn.id.__contains__("ATP_AMP_mc") or rxn.id.__contains__("ATP_ADP_mc") or rxn.id.__contains__("ATP_pc") or rxn.id.__contains__("AMP_ATP_xc") or rxn.id.__contains__("ATP_ADP_Pi_pc"):
        continue
      sto=rxn.metabolites.get(met)
      sto1=rxn.metabolites.get(met1)
      if outfile!="":
        fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution.get(rxn.id)*(sto+sto1))+"\t"+met.compartment+"\n")
      ATPdict[rxn.id]=solution.get(rxn.id)*(sto+sto1)
      if solution.get(rxn.id)*(sto+sto1) > 0:
        total = total + (solution.get(rxn.id)*(sto+sto1))
  if outfile!="":
    fout.close()

  tempDict = dict()
  for rxn in ATPdict.keys():
    tempDict[rxn]=abs(ATPdict[rxn])

  #sort ATPdict by values
  import operator
  sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)

  ATPdict2 = dict()
  ATPdict2["Others-pos"]=0
  ATPdict2["Others-neg"]=0
  baseline = dict()
  pos_base=0
  neg_base=0
  i=0
  for TEMP in sorted_by_value:
    rxn = TEMP[0]
    if ATPdict[rxn]>0:
      if ATPdict[rxn] < total*0.05:
        if percentage:
          ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+float(ATPdict[rxn]*100)/total
        else:
          ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+ATPdict[rxn]
        continue
      base = pos_base
      if percentage:
        ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
        pos_base = pos_base + float(ATPdict[rxn]*100)/total
      else:
        pos_base = pos_base + ATPdict[rxn]
        ATPdict2[rxn]=ATPdict[rxn]
    else:
      if abs(ATPdict[rxn]) < total*0.05:
        if percentage:
          ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+float(ATPdict[rxn]*100)/total
        else:
          ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+ATPdict[rxn]
        continue
      base = neg_base
      if percentage:
        ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
        neg_base = neg_base + float(ATPdict[rxn]*100)/total
      else:
        neg_base = neg_base + ATPdict[rxn]
        ATPdict2[rxn]=ATPdict[rxn]
    i=i+1
    baseline[rxn]=base
  baseline["Others-pos"]=pos_base
  baseline["Others-neg"]=neg_base

  if show_plot:
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    plt.figure(figsize=(3,4))
    for rxn in ATPdict2.keys():
      if colourDict.keys().__contains__(rxn):
        plt.bar(1,ATPdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=colourDict[rxn])
      else:
        plt.bar(1,ATPdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn)
    plt.xlim(0.8,1.2)
    if percentage:
      plt.ylabel("ATP produced/consumed (%)")
    else:
      plt.ylabel("ATP produced/consumed (in moles)")
    handles, labels = plt.gca().get_legend_handles_labels()
    labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
    handles2=[handles[labels.index(i)] for i in labels2]
    lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
    plt.axhline(0,linestyle="--",color="black")
    plt.tight_layout
    plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight')


#####################################################################
# This function generates ATP budgets for a given flux distribution #
# inputs: 1) an FBA model, 2) a dictionary object with reaction ids #
# as keys and reaction fluxes as values, 3) name of output file (op-#
# -tional), 4) Option to show plots, 5) If choosing to show plot, c-#
# -hoose wether to use percentage or absolute values in the plot 6) #
# Provide a day or night indicator tag to specify day or night NAD(-#
# -P)H summary 7) a destination file to save plot to 8) a dictionary#
# to specify colour for fluxes in plot                              #
#####################################################################
def generateNADHNADPHbudget(model,solution,outfile="",show_plot=True,percentage=False,day_or_night_tag="1",save_plot_to="temp",colourDict={}):
    if outfile!="":
        fout = open(outfile,"w")
    Reddict = dict()
    total = 0
    for red in ["NADPH","NADH"]:
        for p in ("c","p","m","x"):
            if len(model.metabolites.query(red+"_"+p+day_or_night_tag))==0:
                continue
            met=model.metabolites.get_by_id(red+"_"+p+day_or_night_tag)
            for rxn in met.reactions:
                sto=rxn.metabolites.get(met)
                sto1=0#rxn.metabolites.get(met1)
                if outfile!="":
                    fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution.get(rxn.id)*(sto+sto1))+"\t"+met.compartment+"\n")
                Reddict[rxn.id]=solution.get(rxn.id)*(sto+sto1)
                if solution.get(rxn.id)*(sto+sto1) > 0:
                    total = total + (solution.get(rxn.id)*(sto+sto1))
    if outfile!="":
        fout.close()

    tempDict = dict()
    for rxn in Reddict.keys():
      tempDict[rxn]=abs(Reddict[rxn])

    #sort by values
    import operator
    sorted_by_value = sorted(tempDict.items(), key= lambda x:x[1],reverse=True)



    Reddict2 = dict()
    Reddict2["Others-pos"]=0
    Reddict2["Others-neg"]=0
    baseline = dict()
    pos_base=0
    neg_base=0
    i=0
    for TEMP in sorted_by_value:
        rxn = TEMP[0]
        if Reddict[rxn]>0:
            if Reddict[rxn] < total*0.05:
                if percentage:
                    Reddict2["Others-pos"]=Reddict2["Others-pos"]+float(Reddict[rxn]*100)/total
                else:
                    Reddict2["Others-pos"]=Reddict2["Others-pos"]+Reddict[rxn]
                continue
            base = pos_base
            if percentage:
                Reddict2[rxn]=float(Reddict[rxn]*100)/total
                pos_base = pos_base + float(Reddict[rxn]*100)/total
            else:
                pos_base = pos_base + Reddict[rxn]
                Reddict2[rxn]=Reddict[rxn]
        else:
            if abs(Reddict[rxn]) < total*0.05:
                if percentage:
                    Reddict2["Others-neg"]=Reddict2["Others-neg"]+float(Reddict[rxn]*100)/total
                else:
                    Reddict2["Others-neg"]=Reddict2["Others-neg"]+Reddict[rxn]
                continue
            base = neg_base
            if percentage:
                Reddict2[rxn]=float(Reddict[rxn]*100)/total
                neg_base = neg_base + float(Reddict[rxn]*100)/total
            else:
                neg_base = neg_base + Reddict[rxn]
                Reddict2[rxn]=Reddict[rxn]
        i=i+1
        baseline[rxn]=base
    baseline["Others-pos"]=pos_base
    baseline["Others-neg"]=neg_base

    if show_plot:
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 10}) #sets a global fontsize
        plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
        plt.rcParams['xtick.major.width'] = 1
        plt.rcParams['ytick.major.size'] = 5
        plt.rcParams['ytick.major.width'] = 1
        plt.rcParams['axes.linewidth']=2 # makes axes line thicker
        plt.figure(figsize=(3,4))
        for rxn in Reddict2.keys():
            if colourDict.keys().__contains__(rxn):
              plt.bar(1,Reddict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn,color=colourDict[rxn])
            else:
              plt.bar(1,Reddict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn)
        plt.xlim(0.8,1.2)
        if percentage:
            plt.ylabel("NAD(P)H produced/consumed (%)")
        else:
            plt.ylabel("NAD(P)H produced/consumed (in moles)")
        handles, labels = plt.gca().get_legend_handles_labels()
        labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
        handles2=[handles[labels.index(i)] for i in labels2]
        lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
        plt.axhline(0,linestyle="--",color="black")
        plt.tight_layout
        plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight')


####################################################################
#This function generates a tab seperated file that can be used with#
#Cytoscape to visualize metabolic flux                             #
#                                                                  #
#inputs: 1) a cobra model with feasible solution 2) the name of the#
#output file  3)  the number of cells in the model (eg: 2 for diel #
#C3 and 4 for diel C4)                                             #
#                                                                  #
####################################################################

def generateFluxMap(cobra_model,solution, outfile,phases = 2):
    import cobra
    #solution = cobra.flux_analysis.parsimonious.optimize_minimal_flux(cobra_model)
    #solution = cobra.flux_analysis.parsimonious.pfba(cobra_model)          #If the previous line returns error comment it out and uncomment this line instead

    #open output file for writing
    f = open(outfile,"w");

    #use rxnSet to identify reaction that have already been processed
    rxnSet = set()

    mult=set()
    #Looping through all reactions in the model
    for rxn in cobra_model.reactions:
        #Get the ID
        RXN=rxn.id
        #declare a boolean variable multFlag to keep track of whether this reaction is present in multiple models
        multFlag=False

        #check if the reaction has already been processed before and if yes skip this run in the loop
        if(rxnSet.__contains__(RXN)):
            continue
        if rxn.id.__contains__("EX") or rxn.id.__contains__("Transfer"):
            multFlag=False
        #check if the reaction ends with one or two i.e it is present more than once in the model
        elif(["1","2","3","4","5","6","7","8","9"].__contains__(rxn.id[len(rxn.id)-1])):
            #change the id to without the suffix 1-9 and declare it as a reaction which has multiple instances
            RXN = rxn.id[0:len(rxn.id)-1]
            multFlag=True
        elif rxn.id[len(rxn.id)-2:] == "10":
            #change the id to without the suffix 10 and declare it as a reaction which has multiple instances
            RXN = rxn.id[0:len(rxn.id)-2]
            multFlag=True

        #if metabolite has multiple instances
        values = dict()
        status1 = dict()
        status2 = dict()
        if(multFlag):
            tempvalue = list()
            temp1 = list()
            temp2 = list()
            mult.add(RXN)
            #add the reaction we are about to process to the reactions processed list
            for i in range(1,phases+1):
                rxnSet.add(RXN+str(i))
                if(round(float(solution.fluxes.get(RXN+str(i)))*10000000) == 0):
                    tempvalue.append(0)
                    temp1.append("none")
                    temp2.append("none")
                elif(float(solution.fluxes.get(RXN+str(i)))*10000 > 0):
                    tempvalue.append(solution.fluxes.get(RXN+str(i)))
                    temp1.append("produced")
                    temp2.append("consumed")
                elif(float(solution.fluxes.get(RXN+str(i)))*10000 < 0):
                    tempvalue.append(solution.fluxes.get(RXN+str(i)))
                    temp1.append("consumed")
                    temp2.append("produced")
            values[RXN] = tempvalue
            status1[RXN] = temp1
            status2[RXN] = temp2

            #select 1 reaction so that we can identify the reactants and products which can be then used to generate the edge shared_name
            rxn=cobra_model.reactions.get_by_id(RXN+"1")

            for reac in rxn.reactants:
                REAC=reac.id
                if(REAC.__contains__("1")):
                    if(REAC.rindex("1")==len(REAC)-1):
                        REAC=REAC[0:len(REAC)-1]
                    f.write("R_"+RXN+" (reaction-reactant) M_"+REAC)
                    for i in range(1,phases+1):
                        f.write("\t"+str(values[RXN][i-1])+"\t"+str(status2[RXN][i-1]))
                    f.write("\n")
                if(RXN.__contains__("biomass")):
                    f.write("R_"+RXN+" (reaction-product)) M_"+REAC)
                    for i in range(1,phases+1):
                        f.write("\t"+str(values[RXN][i-1])+"\t"+str(status1[RXN][i-1]))
                    f.write("\n")
            for prod in rxn.products:
                PROD=prod.id
                if(PROD.__contains__("1")):
                    if(PROD.rindex("1")==len(PROD)-1):
                        PROD=PROD[0:len(PROD)-1]
                f.write("R_"+RXN+" (reaction-product) M_"+PROD)
                for i in range(1,phases+1):
                    f.write("\t"+str(values[RXN][i-1])+"\t"+str(status1[RXN][i-1]))
                f.write("\n")
            if(RXN.__contains__("biomass")):
                f.write("R_"+RXN+" (reaction-reactant) M_"+REAC)
                for i in range(1,phases+1):
                    f.write("\t"+str(values[RXN][i-1])+"\t"+str(status2[RXN][i-1]))
                f.write("\n")
        else:
            #add the reaction we are about to process to the reactions processed list
            rxnSet.add(RXN)
            if(round(float(solution.fluxes.get(rxn.id))*10000000) == 0):
                value = 0;
                status1= "none";
                status0= "none";
            elif(solution.fluxes.get(rxn.id)*10000 > 0):
                value = solution.fluxes.get(rxn.id)*1000;
                status1= "produced";
                status0= "consumed";
            elif(solution.fluxes.get(rxn.id)*10000 < 0):
                value = solution.fluxes.get(rxn.id)*1000;
                status1= "consumed";
                status0= "produced";

            for reac in rxn.reactants:
                REAC=reac.id
                if(REAC.__contains__("1")):
                    if(REAC.rindex("1")==len(REAC)-1):# or (met.id.rindex("2")==len(rxn.id)-1):
                        REAC=REAC[0:len(REAC)-1]
                f.write("R_%s (reaction-reactant) M_%s\t%s\t%s\t0\tnone\n" % (RXN,REAC,value,status0));
            for prod in rxn.products:
                PROD=prod.id
                if(PROD.__contains__("1")):
                    if(PROD.rindex("1")==len(PROD)-1):# or (met.id.rindex("2")==len(rxn.id)-1):
                        PROD=PROD[0:len(PROD)-1]
                f.write("R_%s (reaction-product) M_%s\t%s\t%s\t0\tnone\n" % (RXN,PROD,value,status1));

    f.close();


def writeSolutionFluxesToFile(sol,outfile,model):
    import pandas as pd
    rxnList = list()
    eqnList = list()
    fluxList = list()
    EClist = list()
    for rxn in model.reactions:
        rxnList.append(rxn.id)
        eqnList.append(rxn.reaction)
        if not "PROTEIN CLASS" in rxn.notes.keys():
            EClist.append("")
        else:
            EClist.append(rxn.notes.get("PROTEIN CLASS")[0])
        fluxList.append(sol.fluxes.get(rxn.id))
    df = pd.DataFrame(data={"ID":rxnList,"EC number":EClist,"reaction":eqnList,"flux":fluxList})
    df = df[['ID','EC number', 'reaction', 'flux']]
    df.to_csv(outfile)
    return


def setupMultiphaseModel(model,phases,transferMets):
    from cobra.core import Metabolite, Reaction
    import re

    if phases < 2:
        print("error: number of phases specified is less than 2")
        return

    #if a reaction is A <- B, remake it as B -> A
    for rxn in model.reactions:
        if rxn.lower_bound < 0 and rxn.upper_bound == 0:
            lb = rxn.lower_bound
            for met in rxn.metabolites:
                coeff = rxn.metabolites.get(met)
                rxn.add_metabolites({met:-2*coeff})
                rxn.lower_bound = 0
                rxn.upper_bound = lb


    modelDict = dict()
    #create a copy of all model elements for each phase
    for i in range(1,phases+1):
        temp_model = model.copy()
        for met in temp_model.metabolites:
            met.id = met.id+str(i)
            met.compartment = met.compartment+str(i)
        for rxn in temp_model.reactions:
            rxn.id = rxn.id+str(i)
        modelDict[i] = temp_model
        if i==1:
            final_model = temp_model.copy()
        else:
            final_model = final_model+temp_model
            #add metabolites that might have been missed
            for met in temp_model.metabolites:
                if not final_model.metabolites.__contains__(met.id):
                    final_model.add_metabolites(met.copy())

    #add representation of phloem
    for i in range(1,phases+1):
        rxn = final_model.reactions.get_by_id("Phloem_output_tx"+str(i))
        met = Metabolite("X_Phloem_contribution_t"+str(i),
                         name="X_Phloem_contribution[t]",
                         compartment="t"+str(i))
        rxn.add_metabolites({met:1})

    return final_model


def plotOsmoticContent(cobra_model2,SOL):
    '''
    This function plots osmotic content of the 10-phase developing tomato fruit
    model
    Inputs:1) cobra model 2) solution object
    Outputs:
    '''
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    import numpy as np

    DPA = list(np.linspace(8,59,11))
    for i in range(0,len(DPA)):
        DPA[i]=round(DPA[i],1)

    plt.rcParams.update({'font.size': 15}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth'] = 3 # makes axes line thicker
    plt.figure(figsize=(6,8))

    ax1 = plt.subplot(1,1,1)

    xlist = list()
    ylist = list()
    y2list = list()
    y3list = list()
    y4list = list()
    y5list = list()
    y6list = list()
    y7list = list()
    y8list = list()
    y9list = list()
    y10list = list()
    y11list = list()
    y12list = list()
    y13list = list()
    y14list = list()
    for i in range(1,11):
        xlist.append(DPA[i])
        met = cobra_model2.metabolites.get_by_id("WCO_"+str(i))
        rxn = cobra_model2.reactions.get_by_id("WCOsetter_tx"+str(i))
        #print(rxn.x)
        y1 = rxn.metabolites.get(met)*SOL.fluxes.get(rxn.id)*-1
        ylist.append(y1)
        y2 = 0
        y3 = 0
        y4 = 0
        y5 = 0
        y6 = 0
        y7 = 0
        y8 = 0
        y9 = 0
        y10 = 0
        y11 = 0
        y12 = 0
        y13 = 0
        y14 = 0
        for s in ["VO"]:
            for r in cobra_model2.metabolites.get_by_id(s+"_"+str(i)).reactions:
                if r.metabolites.get(cobra_model2.metabolites.get_by_id(s+"_"+str(i))) > 0:
                    y2 = y2 + (SOL.fluxes.get(r.id)*r.metabolites.get(cobra_model2.metabolites.get_by_id(s+"_"+str(i))))
        y2list.append(y2)
        for s in ["CO"]:
            for r in cobra_model2.metabolites.get_by_id(s+"_"+str(i)).reactions:
                if r.metabolites.get(cobra_model2.metabolites.get_by_id(s+"_"+str(i))) > 0:
                    y3 = y3 + abs(SOL.fluxes.get(r.id)*r.metabolites.get(cobra_model2.metabolites.get_by_id(s+"_"+str(i))))
        y3list.append(y2+y3)
        for s in ["SUCROSE_v","GLC_v","FRU_v"]:
            if i!=10:
                rxn = cobra_model2.reactions.get_by_id(s+"_Transfer"+str(i)+str(i+1))
                y4 = y4 + SOL.fluxes.get(rxn.id)
            else:
                rxn = cobra_model2.reactions.get_by_id("Biomass_"+s+"_tx10")
                y4 = y4 + SOL.fluxes.get(rxn.id)
        y4list.append(y4)
        for s in ["SUCROSE_c","GLC_c","FRU_c"]:
            if i!=10:
                rxn = cobra_model2.reactions.get_by_id(s+"_Transfer"+str(i)+str(i+1))
                y5 = y5 + SOL.fluxes.get(rxn.id)
            else:
                rxn = cobra_model2.reactions.get_by_id("Biomass_"+s+"_tx10")
                y5 = y5 + SOL.fluxes.get(rxn.id)
        y5list.append(y2+y5)
        for s in ["MAL_v","CIT_v"]:
            if i!=10:
                rxn = cobra_model2.reactions.get_by_id(s+"_Transfer"+str(i)+str(i+1))
                y6 = y6 + SOL.fluxes.get(rxn.id)
            else:
                rxn = cobra_model2.reactions.get_by_id("Biomass_"+s+"_tx10")
                y6 = y6 + SOL.fluxes.get(rxn.id)
        y6list.append(y4+y6)
        for s in ["MAL_c","CIT_c"]:
            if i!=10:
                rxn = cobra_model2.reactions.get_by_id(s+"_Transfer"+str(i)+str(i+1))
                y7 = y7 + SOL.fluxes.get(rxn.id)
            else:
                rxn = cobra_model2.reactions.get_by_id("Biomass_"+s+"_tx10")
                y7 = y7 + SOL.fluxes.get(rxn.id)
        y7list.append(y2+y5+y7)
        for s in ["GLN","ASN","SER","GLY","THR","L_ALPHA_ALANINE","4_AMINO_BUTYRATE","VAL","ILE","PHE","LEU","LYS","ARG","L_ASPARTATE","GLT","bHIS","MET","PRO","TRP","TYR","CYS"]:
            #if s == "GLN" or s=="LEU" or s=="CYS":
            #    continue
            if i!=10:
                rxn = cobra_model2.reactions.get_by_id(s+"_v_Transfer"+str(i)+str(i+1))
                y8 = y8 + SOL.fluxes.get(rxn.id)
            else:
                if s=="bHIS":
                    s="HIS"
                rxn = cobra_model2.reactions.get_by_id("Biomass_"+s+"_v_tx10")
                y8 = y8 + SOL.fluxes.get(rxn.id)
        y8list.append(y4+y6+y8)
        for s in ["GLN","ASN","SER","GLY","THR","L_ALPHA_ALANINE","4_AMINO_BUTYRATE","VAL","ILE","PHE","LEU","LYS","ARG","L_ASPARTATE","GLT","HIS","MET","PRO","TRP","TYR","CYS"]:
            #if s == "GLN" or s=="LEU" or s=="CYS":
            #    continue
            if i!=10:
                rxn = cobra_model2.reactions.get_by_id(s+"_c_Transfer"+str(i)+str(i+1))
                y9 = y9 + SOL.fluxes.get(rxn.id)
            else:
                rxn = cobra_model2.reactions.get_by_id("Biomass_"+s+"_c_tx10")
                y9 = y9 + SOL.fluxes.get(rxn.id)
        y9list.append(y2+y5+y7+y9)
        for s in ["STARCH_p"]:
            if i!=10:
                rxn = cobra_model2.reactions.get_by_id(s+"_Transfer"+str(i)+str(i+1))
                y10 = y10 + SOL.fluxes.get(rxn.id)
            else:
                #rxn = cobra_model2.reactions.get_by_id("Biomass_"+s+"_tx10")
                #y10 = y10 + rxn.x
                y10 = 0
        y10list.append(y1+y10)
        for s in ["KI_v","MGII_v","CAII_v","AMMONIUM_v"]:
            if i!=10:
                rxn = cobra_model2.reactions.get_by_id(s+"_Transfer"+str(i)+str(i+1))
                y11 = y11 + SOL.fluxes.get(rxn.id)
            else:
                if s=="KI_v":
                    y11 = y11 + SOL.fluxes.get("K_biomass10")
                elif s=="MGII_v":
                    y11 = y11 + SOL.fluxes.get("Mg_biomass10")
                elif s=="CAII_v":
                    y11 = y11 + SOL.fluxes.get("Ca_biomass10")
                else:
                    y11 = y11 + 0
        y11list.append(y4+y6+y8+y11)
        for s in ["KI","MGII","CAII"]:
            if i!=10:
                y12 = y12 + 0
        #        rxn = cobra_model2.reactions.get_by_id(s+"_Transfer"+str(i)+str(i+1))
        #        y12 = y12 + solution2.fluxes.get(rxn.id)
            else:
                y12 = y12 + SOL.fluxes.get(s+"_biomass_c10")
        y12list.append(y2+y5+y7+y9+y12)
        for s in ["NITRATE_v"]:
            if i!=10:
                rxn = cobra_model2.reactions.get_by_id(s+"_Transfer"+str(i)+str(i+1))
                y13 = y13 + SOL.fluxes.get(rxn.id)
            else:
                y13 = y13 + SOL.fluxes.get("NITRATE_biomass10")
        y13list.append(y4+y6+y8+y11+y13)
        for s in ["NITRATE_c"]:
            if i!=10:
                y14 = y14 + 0
                #rxn = cobra_model2.reactions.get_by_id(s+"_Transfer"+str(i)+str(i+1))
                #y14 = y14 + solution2.fluxes.get(rxn.id)
            else:
                y14 = y14 + SOL.fluxes.get("NITRATE_biomass_c10")
        y14list.append(y2+y5+y7+y9+y12+y14)


    ax1.plot(xlist,y3list,"-",color="black")
    ax1.plot(xlist,y2list,"--",color="black")
    #ax1.fill_between(xlist,y2list,alpha=0.1,color="grey",label="vacuolar contribution",linewidth=0)
    ax1.fill_between(xlist,y4list,alpha=0.3,color="blue",label="vac-SUG",hatch="//",linewidth=0)
    ax1.fill_between(xlist,y4list,y6list,alpha=0.3,color="red",label="vac-OA",hatch="//",linewidth=0)
    ax1.fill_between(xlist,y6list,y8list,alpha=0.3,color="green",label="vac-AA",hatch="//",linewidth=0)
    ax1.fill_between(xlist,y8list,y11list,alpha=0.3,color="purple",label="vac-Cation",hatch="//",linewidth=0)
    ax1.fill_between(xlist,y11list,y13list,alpha=0.3,color="magenta",label="vac-Nitrate",hatch="//",linewidth=0)
    #ax1.fill_between(xlist,y2list,y3list,alpha=0.1,color="white",label="cytosolic contribution",linewidth=0)
    ax1.fill_between(xlist,y2list,y5list,alpha=0.7,color="blue",label="cyt-SUG",hatch="\\\\",linewidth=0)
    ax1.fill_between(xlist,y5list,y7list,alpha=0.7,color="red",label="cyt-OA",hatch="\\\\",linewidth=0)
    ax1.fill_between(xlist,y7list,y9list,alpha=0.7,color="green",label="cyt-AA",hatch="\\\\",linewidth=0)
    ax1.fill_between(xlist,y9list,y12list,alpha=0.7,color="purple",label="cyt-Cation",hatch="\\\\",linewidth=0)
    ax1.fill_between(xlist,y12list,y14list,alpha=0.7,color="magenta",label="cyt-Nitrate",hatch="\\\\",linewidth=0)
    #ax1.fill_between(xlist,y9list,y10list,alpha=0.3,color="black",label="starch",hatch="\\\\",linewidth=0)
    ax1.fill_between(xlist,ylist,y10list,alpha=0.7,color="orange",label="starch (GLC units)",linewidth=0)

    handles, labels = plt.gca().get_legend_handles_labels()
    labels.reverse()
    handles.reverse()
    lgd=plt.legend(handles,labels,bbox_to_anchor=(1.1, 1.02),fontsize=15)
    plt.ylabel("Metabolite content (mmol)")
    plt.xlabel("DPA (days)")
    plt.xlim(DPA[1],DPA[10])
    #plt.xlim(DPA[1],20)
    #plt.ylim(0,1.2)
    plt.show()
