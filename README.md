# Manage kmers and index genomic data with kmtricks, kmindex, and back_to_sequences

Pierre Peterlongo. Ebame-9.

Lecture slides Pierre_Peterlongo_kmers.pdf

## Chapter 0. Preliminaries

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

1. **kmtricks**. As it takes several minutes, you installed kmtricks yesterday. It is accessible by `~/kmtricks/bin/kmtricks`. 

<details><summary>Here is how you installed it:</summary>
<p>
  
```bash
sudo apt-get install cmake --yes
git clone --recursive https://github.com/tlemane/kmtricks
cd kmtricks
./install.sh 
```

- `git clone` needs ~1mn
- `./install.sh` needs ~6mn30

</p>
</details>


Note that kmtricks can also be installed by other means (see this [doc](https://github.com/tlemane/kmtricks/wiki)).

2. **kmindex**. You will install it, see the dedicated section below
  
3. **back_to_sequences**. You will install it, see the dedicated section below
  

## Chapter 1. Testing [kmtricks](https://github.com/tlemane/kmtricks/)

### 1.1 Create a counted kmer matrix

- Create a directory and go there

```bash
mkdir ~/kmtricks_tests && cd ~/kmtricks_tests
```

- Create a file of file containing two datasets: SRR8652861_1 and SRR8653248_1.

```bash
echo "D1: /home/ubuntu/data/public/sra/SRR8652861_1.fastq.gz" > fof_sra.txt
echo "D2: /home/ubuntu/data/public/sra/SRR8653248_1.fastq.gz" >> fof_sra.txt
```

**Question1**: What is the difference between '>' and '>>' in the previous command?

<details><summary>Answer</summary>
<p>

- `cmd > file.txt` enables to store the output of `cmd` in `file.txt`, erasing its potential previous content.
- `cmd >> file.txt` enables to add the output of `cmd` at the end of `file.txt`.
</p>
</details>

**Question2**: What is the cumulated size of `D1`and `D2`?

<details><summary>Answer</summary>
<p>

`ls -lh /home/ubuntu/data/public/sra/SRR8652861_1.fastq.gz /home/ubuntu/data/public/sra/SRR8653248_1.fastq.gz`
This is 255M+185M = 440M

If you want to know the size of uncompressed files: 
- `zcat /home/ubuntu/data/public/sra/SRR8652861_1.fastq.gz | wc -c`
  - 976890579 bytes that are ~932MB (I use https://www.matisse.net/bitcalc for conversions)
- `zcat /home/ubuntu/data/public/sra/SRR8653248_1.fastq.gz | wc -c`
  - 712426733 bytes that are ~679MB.

So the total size of uncompressed files is ~1.5GB

</p>
</details>

**Question3**: How to consider not only reads1, but to consider that `D1`is composed of `SRR8652861_1.fastq.gz` **and** `SRR8652861_2.fastq.gz`, and `D2` is composed of `SRR8653248_1.fastq.gz` **and** `SRR8653248_2.fastq.gz`?

<details><summary>Answer</summary>
<p>
More info here: https://github.com/tlemane/kmtricks/wiki/Input-data
  
With kmindex one line may correspond to multiple input files. So if we wish `D1` to be  `SRR8652861_1.fastq.gz` **and** `SRR8652861_2.fastq.gz` we can indicate this line in the for D1:

`D1: /home/ubuntu/data/public/sra/SRR8652861_1.fastq.gz ; /home/ubuntu/data/public/sra/SRR8652861_2.fastq.gz`

Finally, we can build the fof.txt for D1 and D2 by:

```bash
echo "D1: /home/ubuntu/data/public/sra/SRR8652861_1.fastq.gz ; /home/ubuntu/data/public/sra/SRR8652861_2.fastq.gz" > fof_sra.txt
echo "D2:  /home/ubuntu/data/public/sra/SRR8653248_1.fastq.gz ; /home/ubuntu/data/public/sra/SRR8653248_2.fastq.gz" >> fof_sra.txt
```

(but don't do this, I'd take too much time for this tuto)
</p>
</details>

- Create the raw counting matrix

```bash
~/kmtricks/bin/kmtricks pipeline --file fof_sra.txt --run-dir ./matrix_example --mode kmer:count:bin --hard-min 2 --cpr -t 16
```
(should last ~1mn)

**Question4**: Check and understand all the 7 arguments & options used here:

1. `pipeline`
2. `--file fof_sra.txt`
3. `--run-dir ./matrix_example`
4. `--mode kmer:count:bin`
5. `--hard-min 2`
6. `--cpr`
7. `-t 16`

<details><summary>Answer</summary>
<p>

1. `pipeline`: We used the `pipeline` command from kmtricks. It calls various submodules for you
2. `--file fof_sra.txt`: This is the file of files
3. `--run-dir ./matrix_example`: Where the matrices (among other things) is going to be stored 
4. `--mode kmer:count:bin`: There are many possible output. We chose to represent counted kmers in binary format
5. `--hard-min 2`: We removed kmers having only one occurrence in the input file
6. `--cpr`: Output files are compressed
7. `-t 16`: We use 16 threads
</p>
</details>
  

**Question5**: What is the main output from this command? Can humans read it?
<details><summary>Answer</summary>
<p>

```bash
tree matrix_example/
matrix_example/
|-- build_infos.txt
|-- config_gatb
|   `-- gatb.config
|-- counts
|   |-- partition_0
|   |-- partition_1
|   |-- partition_2
|   `-- partition_3
|-- filters
|-- fpr
|-- hash.info
|-- histograms
|-- howde_index
|-- kmtricks.fof
|-- matrices
|   |-- matrix_0.count.lz4
|   |-- matrix_1.count.lz4
|   |-- matrix_2.count.lz4
|   `-- matrix_3.count.lz4
|-- merge_infos
|   |-- partition0.merge_info
|   |-- partition1.merge_info
|   |-- partition2.merge_info
|   `-- partition3.merge_info
|-- minimizers
|   |-- minimizers.0
|   |-- minimizers.1
|   |-- minimizers.2
|   `-- minimizers.3
|-- options.txt
|-- partition_infos
|   |-- D1.pinfo
|   `-- D2.pinfo
|-- repartition_gatb
|   `-- repartition.minimRepart
|-- run_infos.txt
`-- superkmers
    |-- D1
    |   |-- PartiInfoFile
    |   `-- SuperKmerBinInfoFile
    `-- D2
        |-- PartiInfoFile
        `-- SuperKmerBinInfoFile
```

The most interesting output from this are the four matrices in the `matrix_example/matrices` directory
</p>
</details>

- Dump the matrix in a human-readable format:

```bash
~/kmtricks/bin/kmtricks aggregate --matrix kmer --format text --cpr-in --run-dir matrix_example >  kmer_matrix.txt
```

(should last ~30s)


**Question6**: What is the result and what information does it contain?
<details><summary>Answer</summary>
<p>

Now matrices are readable and can be viewed for instance with `less kmer_matrix.txt`. This file contains kmers having at least 2 occurrences in at least one of the two datasets `D1` and `D2`. 

For instance, the first line `AAAAAAAAAACATTGAACTAATCTAAAAGCA 2 0` indicates that kmer `AAAAAAAAAACATTGAACTAATCTAAAAGCA` has 2 occurrences in `D1` and 0 in `D2`.
</p>
</details>

**Question7**: How many distinct *solid* kmers the two datasets contain?
<details><summary>Answer</summary>
<p>

A _solid_ kmer is a kmer considered as non-erroneous, here having at least 2 occurrences in at least one of the two files.
Just counting the number of lines in the matrix gives the results
```bash
wc -l kmer_matrix.txt
67345405 kmer_matrix.txt
```
So there are ~67 million solid kmers in the two datasets.

</p>
</details>

**Question8 Bonus (If we've time)**: We could have dumped directly the matrices in a human-readable format. How to do this?
<details><summary>Answer</summary>
<p>

kmtricks can output results in a human readable format: 
```bash
~/kmtricks/bin/kmtricks pipeline --file fof_sra.txt --run-dir ./readable_matrix_example --mode kmer:count:text --hard-min 2 --cpr -t 16
```
</p>
</details>

### 1.2 Differences between matrices

Imagine now we are interested in a new dataset SRR8653247_1.fastq.gz. We want kmers neither in SRR8652861_1.fastq.gz nor SRR8653248_1.fastq.gz.

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
<details><summary>Answer</summary>
<p>
  
1. Get the first read: `head -n 1 kmers_only_R1.txt` gives you `AAAAAAAAAACAAAGAGGAGTGGTTTATTAT 2`
2. Check that this kmer `AAAAAAAAAACAAAGAGGAGTGGTTTATTAT` does not occur in `kmer_matrix.txt`: `grep AAAAAAAAAACAAAGAGGAGTGGTTTATTAT  kmer_matrix.txt`
   This gives no results, meaning that the kmer is specific to `R1`
</p>
</details>

**Question10**: How many distinct kmers are only in `R1` and not in `D1` nor `D2`?
<details><summary>Answer</summary>
<p>

Again: `wc -l kmers_only_R1.txt`. There are 19752331 specific to `R1`.
</p>
</details>

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

**Question11**: validate (visually) on the obtained `final_matrix_plugin.txt` that kmers have an abundance of at least 12 and 12 in `D1`and `D2`.
<details><summary>Answer</summary>
<p>

`head final_matrix_plugin.txt` enables to show that first kmers have an abundance at least 12 in the two datasets.
</p>
</details>

**Question12**: Modify the plugin (function `process_kmer`) to retain only kmers seen at least 12 times in `D1` and whose abundance is at least double in `D2`. Recompile, re-test, and validate the output.
<details><summary>Answer</summary>
<p>

In this case we redo exactly the same steps from the start of section 1.4, replacing the `process_kmer` function by
```c++
bool process_kmer(const uint64_t* kmer_data, std::vector<count_type>& count_vector) override
{
// return true if all c[0] > 12 and c[1] > 2*c[0]
  if (count_vector.size() >= 2 && count_vector[0] > m_threshold && count_vector[1] > 2*count_vector[0]){
    return true;
  }
  return false;
}
```

Recompile (still with `-p` option) and re-run kmtricks changing the output directory name: 
```bash
../kmtricks/bin/kmtricks pipeline --plugin ../kmtricks/build/plugins/libmy_plugin.so --plugin-config 12  --file fof_sra.txt --mode kmer:count:bin --hard-min 2 --cpr -t 16 --run-dir matrix_example_plugin2
../kmtricks/bin/kmtricks aggregate --matrix kmer --format text --cpr-in --run-dir matrix_example_plugin2 > final_matrix_plugin2.txt
```

Verify results
```bash
head final_matrix_plugin2.txt 
AAAAAAATATATACGTGAGCAAGAGTTACAT 15 75
AAAAAAATCACCTTGCTTAAGTATAATGTAG 21 50
AAAAAATATATACGTGAGCAAGAGTTACATG 15 74
AAAAAATCACCTTGCTTAAGTATAATGTAGG 21 51
...
```

</p>
</details>

## Chapter 2. Testing [kmindex](https://github.com/tlemane/kmindex)

Up to now, we have created matrices. However, they are used as indexes. This is where kmindex arrives. It can use kmtricks matrices as the basis for building indexes, but it can also handle all the indexing and query pipelines.

### 2.1 Install kmindex

For this lesson, we will use the conda environment. 

```bash
cd
conda create -p kmindex_env --yes
conda activate ./kmindex_env
conda install -c conda-forge -c tlemane kmindex --yes
```
(takes ~1mn30)

Let's go to a specific directory:

```bash
cd
mkdir kmindex_tests && cd kmindex_tests
```

**Question13**: Check that kmindex is installed and have a look at the help.
<details><summary>Answer</summary>
<p>
  
```bash
kmindex -h
```
  Prints the help
</p>
</details>

### 2.2 Create a first index

Let us construct a first index from a set of 100 unitig files from the logan project.


Let's create an index, indexing 100 files each containing unitigs. Files are located in `/ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL`:

First, we need a file of file. We create a file containing, for each .fa file in this directory,: `its accession name: full path`.

```bash
for filename in `ls /ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL/*.fa`; 
  do
  accession=`basename $filename | cut -d "." -f 1`; 
  echo $accession: $filename; 
done > fof.txt
```

**Question14**: Decorticate this for loop, understanding each sub-command.
<details><summary>Answer</summary>
<p>
  
Each `.fa` file is as `/ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL/DRR024015.unitigs.fa`
For each of them, we want to create in the file of file a line as `DRR024015: /ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL/DRR024015.unitigs.fa`

1. The `for` loop enumerates each `.fa` file.
2. `accession='basename $filename | cut -d "." -f 1'`; collects the accession value (`DRR024015` on the example).
  - It removes the path of the filename (removes `/ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL/`)
  - It cuts the remaining (`DRR024015.unitigs.fa`) to keep only what is before first '.' (`cut -d "." -f 1`)
3. it prints the accession followed by ':' and then the full path
4. At the end we store everything in `fof.txt` with `> fof.txt`
</p>
</details>

**Question15**: As files contain unitigs, what should be chosen as the minimum abundance to keep a k-mer? What is the corresponding option name?
<details><summary>Answer</summary>
<p>

All kmers are unique in unitigs. So we should keep all kmers, specifying the threshold to 1.
This is the `--hard-min 1` option
</p>
</details>

**Question16**: The number of distinct kmers per unitig file is below $2^{19}$. We will create bloom filters so that the false positive rate is at most 25%. What should be the size of the bloom filters, knowing that we use 1 hash function?
See for instance [Bloom filter calculator](https://hur.st/bloomfilter/), using a unique hash function, to compute the BF size, knowing the number of kmers to index.


- Note: we used [ntcard](https://github.com/bcgsc/ntCard) to estimate the number of kmers in each file.
<details><summary>Answer</summary>
<p>

https://hur.st/bloomfilter/?n=524288&p=0.25&m=&k=1 

Bloom filter size should be 1822457.


Note: We can also use the following command lines for automatic computation.

```bash
p=0.25
log_bf_size=19
real_bf_size=$( echo "(-( 2^${log_bf_size} )* l($p) / l(2)^2) " | bc -l | cut -d "." -f 1 )
```
</p>
</details>


  

**Question17**: What is the option name used to fix the bloom filter size?
<details><summary>Answer</summary>
<p>

--bloom-size 1822457
</p>
</details>




Now we are ready to run kmindex.

kmindex expects 3 output information:

- the name of the main index (in which several sub-indexes can be registered as we will see later). This is the `--index` parameter.
- the directory name of the constructed index. This is the `--run-dir` parameter.
- the name of the created index. This will be used in the main index. This is the `--register-as` parameter.
  

**Not a Question18**: Run the following command, building an index.


```bash
kmindex build --index index_VRL --run-dir dir_index_VRL --register-as reg_index_VRL --kmer-size 25 --bloom-size 1822457 --hard-min 1 --fof fof.txt
```
(takes 30 seconds)

**Question19**: look at the two created directories, and explain their content.

<details><summary>Answer</summary>
<p>

- One is `dir_index_VRL`. This is where the bloom filter matrices are (in `matrices` sub directory).
- One is `index_VRL`. It contains:
  - a symbolic link to the `dir_index_VRL` directory
  - an `index.json` file that lists the names of indexed files (used later at query time)
</p>
</details>

**Question20**: Use the `index-infos` kmindex command on the created index

<details><summary>Answer</summary>
<p>
```bash
 kmindex index-infos -i index_VRL/ 
```
  Lists the indexed files
</p>
</details>

We arrive now at the most interesting part (?) of the day: let's query our sequences. We want to know where the 3 first reads of `R10000513.unitigs.fa` are eventually among the 100 files we've indexed

Let's first create the query (containing the first 3 unitigs from `ERR10000513unitigs.fa`:

```bash
head -n 6 /ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL/ERR10000513.unitigs.fa > query.fa
```

And now let's run the query:

```bash
kmindex query -i index_VRL -q query.fa -r 0.1 -z 6
```

- Note: `-r 0.1` enables to output indexed datasets for which at least 10% of the kmers from the query are in the dataset.

**Question21**. Was it fast?

<details><summary>Answer</summary>
<p>

YES :)
</p>
</details>

**Question22:** Check and understand the output json file

<details><summary>Answer</summary>
<p>

The previous command created a `output/reg_index_VRL.json` file.
For each queried dataset it shows the ratio of shared kmers between the query and this dataset.

</p>
</details>

**Question23**: verify that each of the three reads was correctly found in the ERR10000513 file

<details><summary>Answer</summary>
<p>

```bash
grep ERR10000513  output/reg_index_VRL.json
```
Enables to show that 100% of kmers of the three unitigs in the query are in set `ERR10000513`.
</p>
</details>

**Question24**: recall the role of the -z parameter, and try different values (including 1 and 8)

<details><summary>Answer</summary>
<p>

We have indexed 25-mers, while our goal was to query 31-mers, so we query 6 consecutive 25-mers per 31-mer. This is the [findere](https://github.com/lrobidou/findere) trick. 
- using `-z 1` the result show more results, that are false positives.
- using `-z 8` we query 33-mers and some results obtained with 31-mers (`-z 6`) show a lower amount of shared k-mers.
</p>
</details>

### 2.3 If we have time: Add a second index to the first one

**Question25** Create a second index with files located in `/ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL2` (create a fof2.txt file of files)
<details><summary>Answer</summary>
<p>

Create a second file of files: 
```bash
for filename in `ls /ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL2/*.fa`;
  do
    accession=`basename $filename | cut -d "." -f 1`;
    echo $accession: $filename;
  done > fof2.txt
```

Now build a second index
```bash
kmindex build --index index_VRL --run-dir dir_index_VRL2 --register-as reg_index_VRL2 --fof fof2.txt --kmer-size 25 --bloom-size 1822457 --hard-min 1 
```
</p>
</details>


**Question26**: understand parameters **index**, **run-dir**, **register-as**

<details><summary>Answer</summary>
<p>

- `index`: this is the main index location
- `run-dir`: this is the run directory where matrices for this dataset (fof2) are going to be stored
- `register-as`: this second index will be _register as_ `reg_index_VRL2` in the main index `index_VRL`
</p>
</details>

**Question27**: understand the structure of the directories

<details><summary>Answer</summary>
<p>

- `dir_index_VRL2` contains the newly created matrices
- `index_VRL` contains now two symbolic links and `index.json` stores the lists the names of the 200 indexed files
</p>
</details>

**Question28**: restart the previous query with the same parameters, did we query both indexes?

<details><summary>Answer</summary>
<p>

```bash
kmindex query -i index_VRL -q query.fa -r 0.1  -z 8
[24/10/22 13:51:27][I:99280] Global index: 'index_VRL'
[24/10/22 13:51:27][I:99280] Starting 'reg_index_VRL' query (100 samples)
[24/10/22 13:51:27][I:99280] Index 'reg_index_VRL' processed. (00s)
[24/10/22 13:51:27][I:99280] Starting 'reg_index_VRL2' query (100 samples)
[24/10/22 13:51:27][I:99280] Index 'reg_index_VRL2' processed. (00s)
[24/10/22 13:51:27][I:99280] Done (00s).
```

We see that the two indexes have been queried
</p>
</details>

### 2.4 [If we have time] Merge a second index to the first one

Even if the two sub-indexes are registered in a unique large one, they are distinct and require two queries. 
With kmindex, we can merge distinct indexes, as long as they were created using the same repartition function.

Let's re-construct a merged index:

First, construct the first one, with distinct names:

```bash
kmindex build --index merged_index_VRL --run-dir dir_index_VRL_merge --register-as reg_index_VRL --kmer-size 25 --bloom-size 1822457 --hard-min 1 --fof fof.txt
```

Second, create the second index, using the repartition computed by the first call:

```bash
kmindex build --index merged_index_VRL --run-dir dir_index_VRL_merge2 --register-as reg_index_VRL2 --from reg_index_VRL -fof fof2.txt --hard-min 1
```
Note that here we reuse all parameters from `merged_index_VRL`, no need to specify bloom size 

Finally, merge the two indexes:

```bash
kmindex merge --index merged_index_VRL --new-name merged_index_VRL --to-merge reg_index_VRL,reg_index_VRL2 --new-path dir_merged_index
```

**Question29**: restart the previous query. How many indexes were queried?

<details><summary>Answer</summary>
<p>

```bash
kmindex query -i merged_index_VRL/ -q query.fa -r 0.1  -z 8
[24/10/22 13:56:16][I:99498] Global index: 'merged_index_VRL'
[24/10/22 13:56:16][I:99498] Starting 'merged_index_VRL' query (199 samples)
[24/10/22 13:56:16][I:99498] Index 'merged_index_VRL' processed. (00s)
[24/10/22 13:56:16][I:99498] Done (00s).
```
We see that we queried here a unique index.

</p>
</details>

## Chapter 3. Find matches between sequences using [back to sequences](https://github.com/pierrepeterlongo/back_to_sequences)

Up to now we "only" found matches between a query and a set of sequences (here unitigs, but it can be genes, reads, ...). Once we have identified a target for a query one may want to go further and find back to which sequences the query is similar.

This is something feasible with [back_to_sequences](https://github.com/pierrepeterlongo/back_to_sequences), a simple but efficient tool.

### 3.1 Install back_to_sequences (b2s)

bs2 is written in Rust. The installation is extremely simple. 
On this brand new VM, on need to install `cargo` (the Rust package manager):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then, source the `cargo` environment: `. "$HOME/.cargo/env"` (or reload the environment).

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

Previously we found a partial match between our query and the unitigs from sample `DRR272392` located in `/ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL/DRR272392.unitigs.fa`.

This was show in output/reg_index_VRL.json, line `DRR272392": 0.4666666666666667,`


We are now going to find which unitigs from `DRR272392` matched the query:

```bash
back_to_sequences --in-kmers query.fa --in-sequences /ifb/data/public/teachdata/ebame/kmindex/SRA_VIRAL/DRR272392.unitigs.fa --out-sequences out_DRR272392_queried.fa --out-kmers out_DRR272392_queried_kmers.txt -m 5
```

**Question30**: What are the two created files, what do they contain?

<details><summary>Answer</summary>
<p>

- file `out_DRR272392_queried.fa` shows the unitigs from `DRR272392.unitigs.fa` that contain at least 5 k-mers from our query.
- file `out_DRR272392_queried_kmers.txt` shows the kmers from our query and their abundance in the contigs (here only 0 or 1 as this are unitigs)
  - in this file we see that only 7 over 40 kmers are in `DRR272392.unitigs.fa`. This is not 46% as output by kmindex. We are here facing some false positives from `kmindex`. In case you don't know what to do this weekend, you may try to re-build the index with a lower FP rate (say 0.1%).
</p>
</details>

**Question31**: Check the possible options of b2s, enabling to obtain more "pseudo-mapping" information.

<details><summary>Answer</summary>
<p>

- adding the `--output-kmer-positions`,the `out_DRR272392_queried_kmers.txt` does not show anymore the number of occurrences of k-mers from the query in the `DRR272392.unitigs.fa`, but instead it shows their occurrences positions in this file (and their orientation)
- adding the `--output-mapping-positions`, the `out_DRR272392_queried.fa` also shows the position of kmers from the query  in unitigs from `DRR272392.unitigs.fa` that contain at least 5 k-mers from our query.
</p>
</details>
