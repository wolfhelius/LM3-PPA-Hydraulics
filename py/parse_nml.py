def process_line(line, intFlag=0):
	temp = line.replace(',','').rstrip().split()
	line = []
	for i in range(5):
		if intFlag:
			line.append(int(temp[i+2]))		
		else:
			line.append(float(temp[i+2]))
	return line

def get_species_meta(nmlfile):

	nml = open(nmlfile, "r")

	species_meta = {}

	for line in nml:    
		if "! species" in line:
			species = process_line(line, 1)
			species_meta["species"] = species
		if "alphaHT" in line:
			alphaHT = process_line(line)
			species_meta["alphaHT"] = alphaHT
		if "thetaHT" in line:
			thetaHT = process_line(line)
			species_meta["thetaHT"] = thetaHT
		if "alphaCA" in line:
			alphaCA = process_line(line)
			species_meta["alphaCA"] = alphaCA
		if "thetaCA" in line:
			thetaCA = process_line(line)
			species_meta["thetaCA"] = thetaCA
		if "alphaCSASW" in line:
			alphaCSASW = process_line(line)
			species_meta["alphaCSASW"] = alphaCSASW
		if "thetaCSASW" in line:
			thetaCSASW = process_line(line)
			species_meta["thetaCSASW"] = thetaCSASW
		if "LMA" in line:
			LMA = process_line(line)
			species_meta["LMA"] = LMA
		if "laimax" in line:
			laimax = process_line(line)
			species_meta["laimax"] = laimax
		if "phiRL" in line:
			phiRL = process_line(line)
			species_meta["phiRL"] = phiRL
		if "maturalage" in line:
			maturalage = process_line(line)
			species_meta["maturalage"] = maturalage
		if "fecundity" in line:
			fecundity = process_line(line)
			species_meta["fecundity"] = fecundity
		if "seedlingsize" in line:
			seedlingsize = process_line(line)
			species_meta["seedlingsize"] = seedlingsize
		if "mortrate_d_c" in line:
			mortrate_d_c = process_line(line)
			species_meta["mortrate_d_c"] = mortrate_d_c
		if "mortrate_d_u" in line:
			mortrate_d_u = process_line(line)
			species_meta["mortrate_d_u"] = mortrate_d_u
		if "rho_wood" in line:
			rho_wood = process_line(line)
			species_meta["rho_wood"] = rho_wood
		if "taperfactor" in line:
			taperfactor = process_line(line)
			species_meta["taperfactor"] = taperfactor
		if "srl" in line:
			srl = process_line(line)
			species_meta["srl"] = srl
		

	return species_meta
