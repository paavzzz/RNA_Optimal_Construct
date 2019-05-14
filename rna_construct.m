# Pavithra Nagarajan
# May 2019
# Purpose of this Project: Generate optimal RNA construct for an RNA of interest (with randomly generated buffer and hairpin regions) that do NOT disrupt native RNA folding.
# MATLAB script that runs my core python script. Required to validate and choose the .py script's generated constructs that have no disruptive folding or base pairing.
# Requires the RNAStructure package from Mathews Lab.
*****This is what must be used TO RUN the rna_construct.py file. ******


#The .py file generates potential RNA constructs that are made up of buffer regions and a hairpin. It stores it in an output file.
#The .m file actually checks to see whether the hairpin and buffer regions do not disrupt the native RNA folding, and also it can be further specified whether the hairpin/buffers should not have ANY small bulges with themselves.
#The final RNA construct sequences are stored in an output file!!

#By default, all the buffers length are 10. The hairpin is length 15.

function rna_construct(rna_of_interest, output_file, hairpin_check, buffer1, buffer2, buffer3, hairpin)
%Examples of what to run in Command Window:
%rna_construct('GAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCA','output.txt','10','1','1','1','1')
%rna_construct('CCCGCCUGGAGGCCGCGGUCGGCCCGGGGCUUCUCCGGAGGCACCCACUGCCACCGCGAAGAGUUGGGCUCUGUCAGCCGCGGG','sample.txt','10','1','1','1','1')

%%%%%%%%
%Example run of Python script on Terminal (in case just want to test python
%file, and not both .py and .m
%rna_of_interest='GAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCA'
%systemCommand = 'python rna_construct.py options.txt GAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCA --hairpincheck=10 --buffer1=1 --buffer2=1 --buffer3=1 --hairpin=1'
%%%%%%%

%This function is required in order to generate the final RNA constructs
%that are ensured to NOT disrupt the native folding of the target RNA of
%interest. User can add specifications on whether the added buffer regions
%must also not base pair, or form part of the hairpin. This is useful
%because sometimes the target RNA has a lot of possibilites in
%buffer/hairpin regions that do not share any partial complements with the
%RNA, and other times it does not. 


%Options for the program.

%output file: REQUIRED
    %for printing potential RNA constructs. This is preliminary, before narrowing down options by running rnastructure on MATLAB. 
    %example: 'options.txt'
    
%rna sequence: REQUIRED
    %example: 'ACUGGGGG'

%hairpincheck: OPTIONAL. default = 10
    %for ensuring hairpin does not base pair with rest of rna, I have put a threshold. I recommend 10.
    %example: '10'. 
    %Helpful Note: I would advise to stay in range of 3 to 10.
    
%buffer1: OPTIONAL. default = 0
    %1 denotes TRUE
    %0 denotes FALSE
    %If you want to restrict all rna constructs to not allow for any base pairing within buffer1, including a small loop or bulge, set to 1.
    %         buffer: XXXXXX  -> dot bracket notation shows base pairing.
    %         outcome: buffer is not used.
    %example:'1'
    %example:'0'
    
%buffer2: OPTIONAL. default = 0
    %refer to buffer1 notes
    %example:'1'
    %example:'0'
    
%buffer3: OPTIONAL. default = 0
    %refer to buffer1 notes
    %example:'1'
    %example:'0'
    
%hairpin: OPTIONAL. default =0
    %1 denotes TRUE
    %0 denoes FALSE
    %allows constructs with only '(((((.....)))))'. No modifications
    %allowed. 
    %           '((((((.....))))))' not allowed (in this case, buffer 1 forms part of the stem in the hairpin.
    %example: '1'
    %example: '0'
    
systemCommand = strcat('python rna_construct.py',{' '}, output_file, {' '},rna_of_interest, {' '},'--hairpincheck=',hairpin_check,{' '},'--buffer1=',buffer1,{' '},'--buffer2=',buffer2,{' '},'--buffer3=',buffer3,{' '},'--hairpin=',hairpin)
[original_structure, bpp_x] = rna_structure(rna_of_interest);
[ans,cmdout]=system(char(systemCommand));
fileID = fopen(output_file); %must match file name specified in systemCommand

rna_constructs = strings;
a = 1;
rna_candidate = fgetl(fileID)
rna_candidate = fgetl(fileID)
cmdout=cmdout(1:4);
b=['sum(structure(1:10)==''.'')==10                                      ';'sum(structure(26:35)==''.'')==10                                     ';'sum(structure((length(structure)-9) : length(structure) )==''.'')==10';'sum(structure(11:25)==''(((((.....)))))'')==15                       ']

while(numel(rna_candidate)>1)

	count = 0
	[structure, bpp_x] = rna_structure(rna_candidate);
	count=0
	debug=1;
    
	for i = [cmdout]
		if str2num(i)==1
			count=count+eval(b(debug,:));

		end
		if str2num(i)==0
			count=count+1;
		end
		debug=debug+1;
	end
	if numel(strfind(structure,original_structure))>0 

		count=count+1;
	end
	if count==5

		rna_constructs(a) = rna_candidate;
    	a = a + 1;
	end
	rna_candidate = fgetl(fileID);
end

fileID = fopen('final_rna_constructs.txt','w');
fprintf(fileID,'%s\n',rna_constructs);
fclose(fileID);
