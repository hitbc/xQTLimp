# impute_tool
## introduction
 Impute_tool is an open source software that implements imputing eQTL information across the genome. The software borrows the method of using the known eQTL statistics to predict unknown eQTL statistics in the ImpG-master software by using the linkage disequilibrium information between SNPs, and completes the imputing of eQTL information. ImpG-master implements the function of completing eQTL imputing in a given window of the user.Based on this, the software reduces the complexity of user input and automatically supplements the information needed to impute the location using 1000G files. In terms of performance, the method of turning on multithreading and preserving LD information is used.
Speed up the impution.
##  Building Impute_tool
 To read gz format, you need zlib.
```bash
sudo apt-get install zlib1g-dev
```
Any C++11 compiler should work.
```bash
make #under impute_tool dir
```
## Usage
### Parameter Description：
-g : the path of Gene_annotation file</br>
-e : the path of Eqtl file</br>
-v : the 100G files dir</br>
-o : the output dir</br>
-t : the num of threads</br>

eg:./Project_gene -g /media/userdisk1/jjpeng/yinquanwei/gencode_v19_gene_annotation.txt -e /media/userdisk1/jjpeng/yinquanwei/Brain_Amygdala.allpairs.txt -v /media/userdisk1/jjpeng/yinquanwei/ -o /media/userdisk1/jjpeng/yinquanwei/output3/ -t 55

