# Manage kmers and index genomic data with kmtricks, kmindex, and back_to_sequences

Pierre Peterlongo. Ebame-9.

Lecture slides Pierre_Peterlongo_kmers.pdf

## 0. Preliminaries

For this tutorial we prepared VMs on wich we installed kmtricks and kmindex. 
After you connect to your VM through SSH run the following command:

```bash
tmux
```

which will start a new session of Tmux (terminal multiplexer).
Tmux is a handy tool for running one or more terminals at the same time, 
with the possibility to detach them (so that they keep running in the background) and/or reattach them to a different terminal.

Detach: `"CTRL-b d"`
Re-attache: `tmux attach`

### Data used for this tutorial.

The data are already on the VM.

1. Data used for testing kmtricks are located in `data/public/sra/`. We will use `SRR8652861_1.fastq.gz`, `SRR8653248_1.fastq.gz`, and `SRR8653247_1.fastq.gz`. They all are Illumina HiSeq 2500 runs of Human infant gut virome.
  
2. Data used for testing kmindex are located in `data/public/teachdata/ebame/kmindex/`, in subdirectories `SRA_VIRAL/` and `SRA_VIRAL2/`. They are **unitigs** from 2x100 various viral RNA datasets. Unitigs were created by the [logan](https://github.com/IndexThePlanet/Logan/blob/main/Unitigs.md) project.
  

### Tools used for this tutorial.

1. **kmtricks**. As it takes several minutes, we already installed kmtricks for you. It is accessible by `~/kmtricks/bin/kmtricks`. If you're interested, here is how it was installed:

```bash
sudo apt-get install cmake --yes
git clone --recursive https://github.com/tlemane/kmtricks
cd kmtricks
./install.sh 
```

Note that kmtricks can also be installed by other means (see this [doc](https://github.com/tlemane/kmtricks/wiki)).

2. **kmindex**. You will install it, see the dedicated section below
  
3. **back_to_sequences**. You will install it, see the dedicated section below
  

## 1. Testing kmtricks

### 1.1 Create a counted kmer matrix

- Create a directory and go there

```bash
mkdir ~/kmtricks_tests && cd kmtricks_tests
```

- Create a file of file containing two datasets: SRR8652861_1 and SRR8653248_1.

```bash
echo "D1: /home/ubuntu/data/public/sra/SRR8652861_1.fastq.gz" > fof_sra.txt
echo "D2: /home/ubuntu/data/public/sra/SRR8653248_1.fastq.gz" >> fof_sra.txt
```

**Question1**: What is the difference between '>' and '>>' in the previous command?

**Question2**: What is the cumulated size of `D1`and `D2`?

**Question3**: How to consider not only reads1, but to consider that `D1`is composed of `SRR8652861_1.fastq.gz` **and** `SRR8652861_2.fastq.gz`, and `D2` is composed of `SRR8653248_1.fastq.gz` **and** `SRR8653248_2.fastq.gz`?

- Create the raw counting matrix

```bash
~/kmtricks/bin/kmtricks pipeline --file fof_sra.txt --run-dir ./matrix_example --mode kmer:count:bin --hard-min 2 --cpr -t 16
```

**Question4**: Check and understand all the 7 arguments & options used here:

1. `pipeline`
2. `--file fof_sra.txt`
3. `--run-dir ./matrix_example`
4. `--mode kmer:count:bin`
5. `--hard-min 2`
6. `--cpr`
7. `-t 16`
  

**Question5**: What is the main output from this command? Can humans read it?

- Dump the matrix in a human-readable format:

```bash
~/kmtricks/bin/kmtricks aggregate --matrix kmer --format text --cpr-in --run-dir matrix_example >  kmer_matrix.txt
```

**Question6**: What is the result and what information does it contain?

**Question7**: How many distinct *solid* kmers the two datasets contain?

**Question8**: We could have dumped directly the matrices in a human-readable format. How to do this?

### 1.2 Differences between matrices

Imagine now we are interested in a new dataset SRR8653247_1.fastq.gz but neither in SRR8652861_1.fastq.gz nor SRR8653248_1.fastq.gz

To do this we can use the filter module. 
However, this module is not installed by default. Let's install it:

```bash
cd ~/kmtricks
./install.sh -m
./bin/kmtricks filter --help
```

Now we can play with this third read set:

```bash
cd ~/kmtricks_tests/
echo "R1: /home/ubuntu/data/public/sra/SRR8653247_1.fastq.gz" > remove_fof.txt
../kmtricks/bin/kmtricks filter --in-matrix matrix_example --key remove_fof.txt --output filtered_matrix_R1_only --hard-min 2 --out-types k --cpr-in
../kmtricks/bin/kmtricks aggregate --run-dir filtered_matrix_R1_only/ --count R1:kmer --format text --output kmers_only_R1.txt
```

**Question9**: validate that first kmer found is only in `R1` and not in `D1` `D2`

**Question10**: How many distinct kmers are only in `R1` and not in `D1` nor `D2`?

### 1.3 Find kmers present in all datasets:

Now let's look at kmers in the three datasets:

```bash
../kmtricks/bin/kmtricks filter --in-matrix matrix_example --key remove_fof.txt --output filtered_matrix_R1_D1_D2 --hard-min 2 --out-types {m,v} --cpr-in
../kmtricks/bin/kmtricks aggregate --matrix kmer --format text  --run-dir filtered_matrix_R1_D1_D2 > kmers_R1_D1_D2.txt
```

### 1.4 Bonus kmtricks: plugins.

Imagine now you're interested in building matrices but only for kmers present at least 12 times in both `D1` and `D2`. In this case, you can filter-out the obtained `kmer_matrix.txt`. However, this can be a heavy file. With kmtricks you can write our plugins, filtering directly the created matrix. Let's have a look.

Create the following plugin in `~/kmtricks/plugins/my_plugin/my_plugin.cpp`

```c++
#include <kmtricks/plugin.hpp>

// DMAX_C is a compile definition set by cmake
using count_type = typename km::selectC<DMAX_C>::type;

class FilterAbundance : public km::IMergePlugin
{
public:
  FilterAbundance() = default;
private:
  unsigned int m_threshold {0};

  // Override process_kmer
  // Discard lines which contain abundances less than a threshold
  bool process_kmer(const uint64_t* kmer_data, std::vector<count_type>& count_vector) override
  {
    for (auto& c : count_vector)
      if (c < m_threshold)
        return false;
    return true;
  }

  // Override configure (not necessary if you don't need configuration)
  // The string is passed to kmtricks with --plugin-config
  // Here it's a simple example where the string is a threshold
  // It could be a path to a config file for instance
  void configure(const std::string& s) override
  {
    m_threshold = std::stoll(s);
  }
};

// Make the plugin loadable
extern "C" std::string plugin_name() { return "filter_abundance"; }
extern "C" int use_template() { return 0; }
extern "C" km::IMergePlugin* create0() { return new FilterAbundance(); }
extern "C" void destroy(km::IMergePlugin* p) { delete p; }
```

The important function here is `process_kmer`. The rest can be forgotten for now.

Let's now compile and run the pipeline, including this plugin:

```bash
cd ~/kmtricks
./install.sh -p
cd ~/kmtricks_tests/
../kmtricks/bin/kmtricks pipeline --plugin ../kmtricks/build/plugins/libmy_plugin.so --plugin-config 12  --file fof_sra.txt --mode kmer:count:bin --hard-min 2 --cpr -t 16 --run-dir matrix_example_plugin
../kmtricks/bin/kmtricks aggregate --matrix kmer --format text --cpr-in --run-dir matrix_example_plugin > final_matrix_plugin.txt
```

**Question11**: validate on the obtained `final_matrix_plugin.txt` that kmers have an abundance of at least 12 and 12 in `D1`and `D2`.

**Question12**: Modify the plugin (function `process_kmer`) to retain only kmers seen at least 12 times in `D1` and whose abundance is at least double in `D2`. Recompile, re-test, and validate the output.

## 2. Testing kmindex

Up to now, we have created matrices. However, they are used as indexes. This is where kmindex arrives. It can use kmtricks matrices as the basis for building indexes, but it can also handle all the indexing and query pipelines.

### 2.1 Install kmindex

For this lesson, we will use the conda environment

```bash
cd
conda create -p kmindex_env --yes
conda activate ./kmindex_env
conda install -c conda-forge -c tlemane kmindex --yes
```

Let's go to a specific directory:

```bash
cd
mkdir kmindex_tests && cd kmindex_tests
```

**Question13**: Check that kmindex is installed and have a look at the help.

### 2.2 Create a first index

Let us construct a first index from a set of 100 unitig files from the logan project.

**Question14**: As files contain unitigs, what should be chosen as the minimum abundance to keep a k-mer? What is the corresponding option name?

**Question15**: The number of distinct kmers per unitig file is below $2^{19}$. We will create bloom filters so that the false positive rate is 25%. What should be the size of the bloom filters?

- Note1: we used [ntcard](https://github.com/bcgsc/ntCard) to estimate the number of kmers in each file.
- Note2: see for instance [Bloom filter calculator](https://hur.st/bloomfilter/), using a unique hash function, to compute the BF size, knowing the number of kmers to index
  

**Question16**: What is the option name used to fix the bloom filter size?

We can also use the following command lines for automatic computation.

```bash
p=0.25
log_bf_size=19
real_bf_size=$( echo "(-( 2^${log_bf_size} )* l($p) / l(2)^2) " | bc -l | cut -d "." -f 1 )
```

Let's create an index, indexing 100 files each containing unitigs. Files are located in `/ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL`:

First, we need a file of file. We create a file containing, for each .fa file in this directory,: `its accession name: full path`.

```bash
for filename in `ls /ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL/*.fa`; 
  do accession=`basename $filename | cut -d "." -f 1`; 
  echo $accession: $filename; 
done > fof.txt
```

**Question17**: Decorticate this for loop, understanding each sub-command.

Now we are ready to run kmindex.

kmindex expects 3 output information:

- the name of the main index (in which several sub-indexes can be registered as we will see later). This is the `--index` parameter.
- the directory name of the constructed index. This is the `--run-dir` parameter.
- the name of the created index. This will be used in the main index. This is the `--register-as` parameter.
  

**Question18**: add the missing options and run the completed command:

```bash
kmindex build --index index_VRL --run-dir dir_index_VRL --register-as reg_index_VRL --kmer-size 25 XXX
```

**Question19**: look at the two created directories, and explain their content.

**Question20**: Use the `index-infos` kmindex command on the created index

We arrive now at the most interesting part (?) of the day: let's query our sequences. We want to know where the 3 first reads of `R10000513.unitigs.fa` are eventually among the 100 files we've indexed

Let's first create the query:

```bash
head -n 6 /ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL/ERR10000513.unitigs.fa > query.fa
```

And now let's run the query:

```bash
kmindex query -i index_VRL -q query.fa -z 5
```

**Question21**. Was it fast?

**Question22:** Check and understand the output json file

**Question23**: verify that each of the three reads was correctly found in the ERR10000513 file

**Question24**: recall the role of the -z parameter, and try different values (including 1 and 8)

### 2.3 Add a second index to the first one

**Question25** Create a second index with files located in `~/data/mydatalocal/SRA_VIRAL2` (create a fof2.txt file of files)

```bash
kmindex build --index index_VRL --run-dir dir_index_VRL2 --register-as reg_index_VRL2 --fof fof2.txt --kmer-size 25 XXX
```

**Question26**: understand parameters **index**, **run-dir**, **register-as**

**Question27**: understand the structure of the directories

**Question28**: restart the previous query with the same parameters, did we query both indexes?

### 2.4 Merge a second index to the first one

Even if the two sub-indexes are registered in a unique large one, they are distinct and require two queries. 
With kmindex, we can merge distinct indexes, as long as they were created using the same repartition function.

Let's re-construct a merged index:

First, construct the first one, with distinct names:

```bash
kmindex build --index merged_index_VRL --run-dir dir_index_VRL_merge --register-as reg_index_VRL --kmer-size 25 XXX
```

Second, create the second index, using the repartition computed by the first call:

```bash
kmindex build --index merged_index_VRL --run-dir dir_index_VRL_merge2 --register-as reg_index_VRL2 --kmer-size 25 --from reg_index_VRL XXX
```

Finally, merge the two indexes:

```bash
kmindex merge --index merged_index_VRL --new-name merged_index_VRL --to-merge reg_index_VRL,reg_index_VRL2 --new-path dir_merged_index
```

**Question29**: restart the previous query. How many indexes were queried?

## 3. Find matches between sequences.

Up to now we "only" found matches between a query and a set of sequences (here unitigs, but it can be genes, reads, ...). Once we have identified a target for a query one may want to go further and find back to which sequences the query is similar.

This is something feasible with [back_to_sequences](https://github.com/pierrepeterlongo/back_to_sequences), a simple but efficient tool.

### 3.1 Install back_to_sequences (b2s)

bs2 is written in Rust. The installation is extremely simple. 
On this brand new VM, on need to install `cargo` (the Rust package manager):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then, source the `cargo` environment: `. "$HOME/.cargo/env"` (or reload the shell).

Once `cargo` is installed, it can be used to install b2s and all its dependencies:

```bash
git clone https://github.com/pierrepeterlongo/back_to_sequences.git
cd back_to_sequences
RUSTFLAGS="-C target-cpu=native" cargo install --path .
cd ..
rm -rf back_to_sequences
```

Now that b2s is installed, just check the help: ` back_to_sequences --help`

### 3.2 Run back_to_sequences

Previously we found a partial match between our query and the unitigs from sample `ERR10000176` located in `/ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL/ERR10000176.unitigs.fa`.

We are now going to find which unitigs from `ERR10000176` matched the query:

```bash
back_to_sequences --in-kmers query.fa --in-sequences /ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL/ERR10000176.unitigs.fa --out-sequences out_ERR10000176_queried.fa --out-kmers out_ERR10000176_queried_kmers.txt -m 5
```

**Question30**: What are the two created files, what do they contain?

**Question31**: Check the possible options of b2s, enabling to obtain more "pseudo-mapping" informations.
