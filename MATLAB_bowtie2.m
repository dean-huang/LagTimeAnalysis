%% Using MATLAB to run bowtie2

%% Things to do before running analysis.

% Download data:
% Download the FASTQ file. Although the data is technically from paired end
% reads, there is no need to use the paired ends. 

% Decompressing data (skip if data isn't zipped):
% First extract the .zip file into your data files. Then use the terminal
% to cd into the parent directory of all the data folders. Then use the
% following command to decompress all the .gz files within the folders: 
% gzip -dv ./*/*.gz
% The flag -d tells it to decompress. The flag -v makes it verbose so that
% you can see which files are being decompressed. 

% Change directory to main workspace.
dir_main = [fileparts(which('MATLAB_bowtie2.m')),'/'];
cd(dir_main);

% Set up directories for the reference files.
ref_strain = 'MG1655';
dir_ref = [dir_main,'ReferenceGenomes/',ref_strain,'/'];
ref_fasta = [ref_strain,'.fasta'];
ref_index = [ref_strain];

% Set up directories for the files.
dir_study = 'Data/Ecoli/LB_Exp/';
filename = 'ERR2403099';

% Name the SAM file output.
samfile_base = '3099';
samfile = ['sam_',samfile_base,'.sam'];


%% Building the reference genome index using bowtie2build. Can skip if already made.

% Makes an index file of the reference genome (in this case MG1655) using
% bowtie2build. This will then be used to map the data onto the reference.
% If it returns 0, that means the build was successful.
cd(dir_ref);
bowtie2build(ref_fasta,ref_index)
cd(dir_main);


%% Mapping a sequence using bowtie2

tic

% Get directory name of the FASTQ file as char data structure.
dir_fastq = [dir_main,dir_study];
cd(dir_fastq);
reads = [dir_fastq,dir([filename,'*.fastq']).name];

% Set extra align options for bowtie2.
alignOpt = Bowtie2AlignOptions;
% bowtie2 will find the best alignment. If you want it to find all
% alignments and sort them by alignment score, do 'All'.
alignOpt.NumAlignments = 'Best';
% The samfile will exclude unaligned reads
alignOpt.ExcludeUnaligned = true;

% Unpaired reads. Set reads2 as empty "". Takes about 7 mins for 2 GB FASTQ
% file. The alignOpt flag includes alignment options specified above.
bowtie2([dir_ref,ref_index],reads,"",samfile,alignOpt);

toc


%% Get SAM file information. Uses samtools.

cd(dir_fastq);

% Set the variable samfile if not already specified.
if ~exist('samfile','var')
    samfile_base = input('Please input the SAM file indicator name with single quotes: ');
    samfile = ['sam_',samfile_base,'.sam'];
    disp('SAM file updated.')
end

% Set the number of bins that the histogram will have. Create a hist_counts
% object that keeps track of the number of reads in each bin for successive
% blocks.
L = 4641652;
nbin = 1000;
edges = linspace(0,L/2,nbin+1);

tic


disp(['Analyzing SAM file for ',samfile_base])
    
system(['samtools view ',samfile,' | cut -f4 > pos_',samfile_base,'.txt']);
pos = load(['pos_',samfile_base,'.txt']);
save(['pos_',samfile_base,'.mat'], 'pos');
delete(['pos_',samfile_base,'.txt']);

toc

% Save the histogram.
hist_counts = histcounts(nonzeros(pos),edges);
plot(1:nbin,hist_counts,'.','MarkerSize',3)
ax = gca;
savefig(['hist_',samfile_base,'.fig'])
exportgraphics(ax,['hist_',samfile_base,'.eps'],'ContentType','vector')


%% Save files 

% Save the histogram counts for future use.
save(['histcounts_',samfile_base,'.mat'],'hist_counts')

% Save the plot as an .eps file for future use.
ax = gca;
exportgraphics(ax,['hist_',samfile_base,'.eps'],'ContentType','vector')
savefig(['hist_',samfile_base,'.fig'])


