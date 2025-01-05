







def best_result(file_name, monomer, rec_lig, receptor, ligand):

	# Function to generate the "best docking results", being the result with the best score and with the residue with the best contact frequency

	file_name_dir = str('./results_2/' + receptor + '_' + ligand + '_folder_'  + str(date.today()) + '/' + receptor + '_' + monomer + '_' + ligand + '/result/') #directory for the docking results
	file_name_path = str(file_name_dir + file_name[:-20] + '.pdb') #directory for the result, identifies as the best result
	des1 = file_name_dir + 'best_docking_results_for_' + file_name[:-24] + '.pdb' #destination directory for the best_docking_result file
	shutil.copyfile(file_name_path,des1)

	#Same thing done with the ligand file only
	ori2='./results_2/' + receptor + '_' + ligand + '_folder_' + str(date.today()) + '/' + receptor + '_' + monomer + '_' + ligand + '/ligand_reserved_pdb/' + file_name
	des2='./results_2/' + receptor + '_' + ligand + '_folder_' + str(date.today()) + '/' + receptor + '_' + monomer + '_' + ligand + '/ligand_reserved_pdb/best_docking_results.pdb'
	shutil.copyfile(ori2,des2)


	# This is to create a copy of that file with 'Z' as the name of the chain in the ligand,
	# it is important for the 3dsjmol visualization
	# (I believe... it was in Yuexin original script and if i remove this part the 3d visualization in the final html has issues)

	with open(str(file_name_dir + 'best_docking_results_for_' + file_name[:-24] + '.pdb'), 'r') as file: #to not modify the chain name for the protein chains
		lines = file.readlines()
		subpart1 = lines[:lines.index(
			'REMARK    Docked ligand coordinates...\n')]  #subpart 1 is from start to 1st line in ligand coordinates
		subpart2 = lines[lines.index(
			'REMARK    Docked ligand coordinates...\n'):] #subpart 2 from 1st line in ligand coordinates to end of file
	with open(str(file_name_dir + 'best_docking_results_for_' + receptor + '_' + monomer + '_' + ligand + '.pdb'), 'w') as file:
		for l in subpart1:
			file.write(l)
		for line in subpart2:
			if line[0:4] == 'ATOM' or line[:6] == 'HETATM' or line[:3] == 'TER':
				newline = line[:21] + 'Z' + line[22:]
				file.write(newline)
			else:
				file.write(line)
	print('best docking result file is generated for ' + file_name[:-24])



def separate_results(monomer, file_dir, first_file_name, dir_final, monomers_list):

	# Function to separate the multimer file into its monomers for every result file created by hex

	ends = [] #this list will be modified with the indices of every monomer's terminal line + the first coordinate's line index
	# Open the .pdb file to separate
	with open (file_dir + first_file_name, 'r+') as r:
		lines = r.readlines()
		for l in lines:
			if l.startswith('ATOM      1  '):
				ends.append(lines.index(l)) #and save the index of the first coordinate's line in the list ends

		# Searches the .pdb files for the lines that indicate the end of a chain
		for l in lines:
			if l[0:3] == 'TER':
				ends.append(lines.index(l)) #and add their indexes in the ends list

		if os.path.isdir(dir_final) == False: #create folder to dump the new monomer file or files
			os.makedirs(dir_final)

		# LOGIC:The end of the previous chain is the start of the current one,
		start_pos = ends[monomers_list.index(monomer)]
		end_pos = ends[monomers_list.index(monomer)+1]

	# or that it is in the range of the monomer we want to isolate
	# It copies every line that is not referencing an atom coordinates
	file_list = os.listdir(file_dir)
	for r in file_list: #for every result file:
		file_path = str(file_dir + '/' + r)
		new_file_path = str(dir_final + r[:-4] + '_' + monomer + '.pdb') #create a new result file which will include only one protein chain, not all
		with open(file_path, 'r') as file:
			lines = [line for line in file.readlines()]
			# Dump in the new file everything before the first coordinate line + between the lines that contain
			# the monomer coordinates + after the last receptor's coordinates
			lines = lines[:ends[1]] + lines[start_pos:end_pos] + lines[ends[-1]:]
		with open(new_file_path, 'w') as file:
			file.writelines(lines)




def separate_monomers(monomer, file_dir, file_name, dir_final, monomers_list):

	# Function to separate the original protein pdb file in its monomers

	# Open the .pdb file to separate
	with open (file_dir + '/' + file_name + '.pdb', 'r+') as r:
		lines = r.readlines()
		ends = [0]

		# Searches the .pdb files for the lines that indicate the end of a chain
		for l in lines:
			if l[0:3] == 'TER':
				ends.append(lines.index(l))
		if os.path.isdir(dir_final) == False:
			os.makedirs(dir_final)
		monomer_pdb = open(dir_final + '/' + file_name + '_' + monomer + '.pdb', 'a+')


		# The end of the previous chain is the start of the current one,
		# 0 was previously included in the list ends to be the start of the first chain
		start_pos = ends[monomers_list.index(monomer)]
		end_pos = ends[monomers_list.index(monomer)+1]

		# It copies every line that is not referencing an atom coordinates
		# or that it is in the range of the monomer we want to isolate
		for l in lines:
			if l[0:4] != 'ATOM' or lines.index(l) in range(start_pos, end_pos):
				monomer_pdb.write(l)
			# It needs to copy also the ligand data (if there is any) which is labeled with SDF
			elif l[17:20] == 'SDF':
				monomer_pdb.write(l)





def ligand_reserved(monomer, rec_lig, receptor, ligand):

	# Function to separate the ligand coordinates of every solution, it's useful to simply the calculation of the contact frequencies

	dir_path = str('./results_2/'+ rec_lig + '_folder_' + str(date.today()) + '/' + receptor + '_' + monomer + '_' + ligand + '/result') #results directory
	print('Isolating ' + rec_lig + '_' + monomer)

	os.makedirs('./results_2/'+ rec_lig + '_folder_' + str(date.today()) + '/' + receptor + '_' + monomer + '_' + ligand + '/ligand_reserved_pdb') #ligand_reserved directory
	file_list = os.listdir(dir_path)
	result_list = []

	# Some operative system will create hidden files, the script consider .pdb files only
	for i in file_list:
		if i[0] != '.' and len(i.split('.')) == 2 and i.split('.')[1] == 'pdb':
			result_list.append(i)
	for r in result_list:
		file_path = str(dir_path + '/' + r)
		ligand_reserved_file_path = str('./results_2/'+ rec_lig + '_folder_' + str(date.today()) + '/' + receptor + '_' + monomer + '_' + ligand + '/ligand_reserved_pdb/' + r[:-4] + '_ligand_reserved.pdb')
		with open(file_path, 'r') as file:
			lines = [line for line in file.readlines()]
			# Everything below the line 'REMARK    Docked ligand coordinates...' is data of the ligand
			lines = lines[lines.index('REMARK    Docked ligand coordinates...\n'):]
		with open(ligand_reserved_file_path, 'w') as file:
			file.writelines(lines)





def result_dict_generator(threshold, monomer, rec_lig, receptor, ligand):

	# Function to calculate the contact frequencies of every amino acid

	result_dir_path = str('./results_2/'+ rec_lig + '_folder_' + str(date.today()) + '/' + receptor + '_' + monomer + '_' + ligand + '/ligand_reserved_pdb/') #directory for the results files, the ligand only ones we created with the ligand_reserved function!
	receptor_file_path = str('./results/receptor_to_dock/monomers_' + receptor + '/' + receptor + '_' + monomer + '.pdb') #directory for the receptor protein pdb file

	# Store every receptor's atom coordinates information as a nested dictionary called 'reference'
	with open(receptor_file_path, 'r') as file:
		reference = {}
		for line in file.readlines():
			if line[0:4] == 'ATOM':
				if int(line[22:27]) in reference:
					reference[int(line[22:27])][int(line[6:11])] = tuple(map(float, filter(None, line[31:54].split(' '))))
				else:
					reference[int(line[22:27])] = {int(line[6:11]) : tuple(map(float, filter(None, line[31:54].split(' '))))}

	#so the reference is {residue: {atom :(x, y, z)}}

	# The energy for each reference element will be stored in dictionary 'ac'
	ac = {}
	file_list = os.listdir(result_dir_path)
	result_list = []

	# Generate the list for all .pdb names in the directory
	for i in file_list:
		if i[0] != '.' and len(i.split('.')) == 2 and i.split('.')[1] == 'pdb':
			result_list.append(i)

	en_list = [] #future list of energies
	file_names = [] #future list of file names
	resi_list = [] #future list of aa

	#reading the first file and saving its lines will make things much quicker for the rest of them
	first_file_path = str(result_dir_path + receptor + '_' + ligand + '0001_' + monomer + '_ligand_reserved.pdb')
	z=open(first_file_path)
	lines_first=z.readlines()
	x=lines_first[2]
	print (x)


	# Store energy values for each ligand_reserved file
	for r in result_list:
		print('current file:' + r)
		energy = ''
		file_path = str(result_dir_path + r)

		with open(file_path) as file:
			lines = file.readlines()
			for l in lines:
				if 'REMARK' in l.split(' ') and 'Energy' in l.split(' '):
					# The energy is divided by the number of results to
					# later obtain an average energy when we will sum the
					energy = (float(l.split(' ')[6][:-1]))/(len(result_list))
					# Generate file and energy list by order
					file_names.append(str(r))
					en_list.append(energy)

			# Go over every coordinate of atoms in the ligand_reserved file and store into coor
			coor = [tuple(map(float, filter(None, line[31:54].split(' '))))
					for line in lines if line[0:4] == 'ATOM']
			lst = []

			for res in reference.keys(): # for each amino acid in the receptor file:
				distances = []

				for atom in coor: # for each atom of the ligand

					for aa in reference[res].keys(): # for each atom of that amino acid
						# check if the distance between atoms of the ligands
						# and of the amino acid are lower than chosen threshold (5)
						distances.append(math.sqrt((reference[res][aa][0] - atom[0]) ** 2 + (reference[res][aa][1] - atom[1])** 2
									 + (reference[res][aa][2] - atom[2]) ** 2))

				if all(d >= threshold for d in distances): #if none of the distances is lower than the threshold, skip
					continue

				else: # if at least one distance is lower then add this aminoacid to the ac dict
					if res in ac.keys():
						ac[res] += energy	# adding energy (previosly divided by the number of results) more times if
					else:				 	# found multiple times, that way you would have an average
						ac[res] = energy

					# Store the resi number into lst
				if res not in lst:
						lst.append(res)
			# Store rei_num for one file into resi_list as a list
			resi_list.append(lst)



	best_result_name = ''
	# Find the resi number with the lowest energy
	red_resi = ''
	for k, v in ac.items():
		if v == min(ac.values()):
			red_resi = k
	print('best_residue: ' + str(red_resi))

	# Find the file that both satisfies the lowest energy and containing the lowest energy resi
	max_en = 0
	for f in file_names:
		if en_list[file_names.index(f)] <= max_en:
			temp = resi_list[file_names.index(f)]
			for i in temp:
				if i == red_resi:
					best_result_name = f


	res_dict_path = result_dir_path + 'res_dict.json'

	# Use the result file from /result/, change the name to best docking result, and convert it into chain Z
	try:
		best_result(best_result_name, monomer, rec_lig, receptor, ligand)
	# sometimes the simulations results are not good enough to satisfy both requirements,
	# it's common especially when one monomer is never close to the ligand.
	# Not including this line would stop an otherwise useful simulation
	except FileNotFoundError:
		f_file = receptor + '_' + ligand + '0001_' + monomer + '_ligand_reserved.pdb'
		best_result(f_file, monomer, rec_lig, receptor, ligand)

	print(ac)

	with open(res_dict_path, 'w') as file:
		file.write(json.dumps(ac))
	print('res_dict.json is generated')
	return ac



#

def color_surfaces(monomer, receptor, ligand, rec_lig):

	# Function to create the nested dictionary with every monomer as key with value a dictionary with its amino acids as keys and contact frequencies as values

	result_dict = {} #this will be the dictionary

	folder_name = str(receptor + '_' + monomer + '_' + ligand)

	if receptor + '_' + monomer not in result_dict.keys():
		result_dict[receptor + '_' + monomer] = {}
	if os.path.isfile('./results_2/' + rec_lig + '_folder_' + str(date.today()) + '/' + folder_name + '/ligand_reserved_pdb/res_dict.json') == False:
		result_dict[receptor+ '_' + monomer][ligand] = result_dict_generator(5, monomer, rec_lig, receptor, ligand)
	else:
		result_dict[receptor+ '_' + monomer][ligand] = eval(
			open('./results_2/' + rec_lig + '_folder_' + str(date.today()) + '/' + folder_name + '/ligand_reserved_pdb/res_dict.json', 'r').read())
		print('res_dict.json previously exists and has read')

	resultjson_path = './results_2/' + rec_lig + '_folder_' + str(date.today()) + '/' + folder_name + '/results.json'

	# Initialize results.json
	ini = {}
	with open(resultjson_path, 'w') as file:
		file.write(json.dumps(ini))
	results = {}
	for r in result_dict: #result_dict is where we have our contact freuquencies
		if r in results.keys():
			for v in result_dict[r]:
				results[r][v] == result_dict[r][v]
		else:
			results[r] = result_dict[r]
	with open(resultjson_path, 'w') as file:
		file.write(json.dumps(results))
	print('result.json is finished')








def pipeline(rec_lig, is_monomer, receptor, ligand, monomers_list):

	print('Current pair:' + rec_lig)

	today_dir = './results_2/' + rec_lig + '_folder_' + str(date.today()) + '/'

	datetoday = str(date.today())

	results_dir = today_dir + rec_lig + '/result/'
	os.makedirs(results_dir)

	#hex_docking(rec_lig, rec_lig, receptor, ligand, datetoday) # CALL HEX

	#results_list = os.listdir(results_dir)
	#first_file_name = str(receptor + '_' + ligand + '0001.pdb')


	# Repeats the analysis for every monomer in the receptor file
	for monomer in monomers_list:
		dir_final = today_dir + receptor + '_' + monomer + '_' + ligand + '/result/'
		print('plotting monomer: ' + monomer + ' with the ligand: ' + ligand)
		separate_results(monomer, results_dir, first_file_name, dir_final, monomers_list)
		ligand_reserved(monomer, rec_lig, receptor, ligand)
		print('Ligands are now reserved in docking results.')
		color_surfaces(monomer, receptor, ligand, rec_lig)
		#plot_frequencies(monomer)





def start():

	# Check if the receptor is a monomer or a complex and save the receptor and ligand names as variables

	receptor_folder = '/home/metyumelkonyan/BCB330/results/receptor_to_dock/'
	receptor_folder_list = os.listdir(receptor_folder)
	ligand_folder = os.listdir('/home/metyumelkonyan/BCB330/results/ligand_to_dock/')
	receptor = "bri1.pdb"
	ligand = "brass.pdb"
	# To check if the receptor is a monomer or not, the script will search the .pdb file
	# for the line that indicated the presence of multiple chains,
	with open(receptor_folder + '/' + receptor, 'r+') as f:
		is_monomer = True
		for x in f.readlines():
			if re.match(r'COMPND   \d CHAIN: \w, \w*', x) != None:
				is_monomer = False
				#if the receptor would be a monomer the regex would be r'COMPND   \d CHAIN: \w;'

				# To make a list of the monomers' labels
				print(receptor + ' identified as a protein complex')
				if x[11:16] == 'CHAIN':
					monomers_list = x.split(': ')[-1].split(', ')
					# The COMPND line ends with ';' therefore it needs to be removed from the last label
					monomers_list[-1] = monomers_list[-1][0]

	rec_lig = receptor + '_' + ligand
	home='/home/metyumelkonyan/BCB330/'
	# To save the terminal output later (very important)
	stdoutOrigin=sys.stdout
	sys.stdout = open(home + 'results_2/Terminal_recordings/' + receptor + '_' + ligand + '_' + str(date.today()) + '.txt' , "w")

	# Call to the pipeline with different parameters whether the receptor is a monomer or a complex
	if is_monomer == False:
		dir_final = './results/receptor_to_dock/monomers_' + receptor
		for monomer in monomers_list:
			print('separating monomer: ' + monomer)
			separate_monomers(monomer, receptor_folder, receptor, dir_final, monomers_list) # To separate the monomers in the multimer file

		pipeline(rec_lig, is_monomer, receptor, ligand, monomers_list)
	else:
		dir_final = './results/receptor_to_dock/monomers_' + receptor
		monomers_list = ['monomer']
		separate_monomers('monomer', receptor_folder, receptor, dir_final, monomers_list) # To analyze the data from hex you still need to separate it.
																						  # It allows to use the same functions in both cases
		pipeline(rec_lig, is_monomer, receptor, ligand, monomers_list)

	#to put together the json files with all the data from all monomers
	new_json = './results_2/'+ rec_lig + '_folder_'  + str(date.today()) + '/' + '/final.json'
	final_json = {}
	min_values = []
	max_values = []

	import json
from datetime import date

# Assuming `monomers_list`, `rec_lig`, `receptor`, `ligand`, `min_values`, and `max_values` are defined elsewhere
for monomer in monomers_list:
    monomer_json = './results_2/' + rec_lig + '_folder_' + str(date.today()) + '/' + str(receptor + '_' + monomer + '_' + ligand) + '/results.json'
    with open(monomer_json, 'r') as file:
        monomer_dict = json.load(file)
        print(monomer_dict)

        monomer_key = list(monomer_dict.keys())[0]
        ligand_key = list(monomer_dict[monomer_key].keys())[0]
        inside_dict = monomer_dict[monomer_key][ligand_key]

        if monomer_dict[monomer_key][ligand_key] == {}:
            continue
        else:
            mini = min(inside_dict.values())
            maxi = max(inside_dict.values())

            min_values.append(mini)
            max_values.append(maxi)



print(min_values)
print(max_values)

abs_max = max(max_values)
abs_min = min(min_values)

print(abs_min)
print(abs_max)


monomer_json = './results_2/' + rec_lig + '_folder_' + str(date.today()) + '/' + str(receptor + '_' + monomer + '_' + ligand) + '/results.json'
with open(monomer_json, 'r') as file:
		monomer_dict = json.load(file)

		monomer_key = list(monomer_dict.keys())[0]

		ligand_key = list(monomer_dict[monomer_key].keys())[0]

		inside_dict = monomer_dict[monomer_key][ligand_key]

		normalized_mon_dict = {monomer_key: {ligand_key:{k: (v-abs_min)/(abs_max - abs_min) for k, v in inside_dict.items()}}}

		final_json.update(normalized_mon_dict)

with open(new_json, 'w') as  file:
		file.write(json.dumps(final_json))
print('result.json is finished')

def	best_normalization(monomer_json,receptor,ligand,final_json):


	# To remove the temporary files and directories created

	shutil.rmtree('./results/receptor_to_dock/monomers_' + receptor, ignore_errors = True)
	shutil.rmtree('./results/receptor_to_dock/' + receptor + 'pdb', ignore_errors = True)
	shutil.rmtree('./results/ligand_to_dock/' + ligand + 'sdf', ignore_errors = True)

	sys.stdout.close()
	sys.stdout=stdoutOrigin
