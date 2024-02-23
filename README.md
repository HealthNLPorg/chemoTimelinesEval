# Evaluation code for ChemoTimelines 2024 

##  Evaluation script versioning

#### Current timestamp (YYYY-MM-DD): 2024-02-23
Update for 2024-02-23 version: Fix `divided by zero exception` Line 405 at `eval_timeline.py`. 
<br>Update for 2024-01-16 version: Initial release.

#### Notice:
We have internally reviewed this script multiple times. However, should there be any concerns or feedback (i.e., if you find errors in this code), please let us know. We are open to feedback on the code until ~February 16, 2024~ **March 10, 2024**, and will respond within 3 calendar days.

~February 21, 2024~ **March 20, 2024** (Updated to reflect the changes in Important Dates), is the cut-off date for the final update in the unlikely event of updating this script. Please mark this date for the final sync (`git pull origin current`) of this repository. The code will not be updated after this date to minimize confusion.


## Summarizing Docker output

Obtain a summarized JSON timeline via running `docker_output_to_timeline.py` on the Docker output TSV via:
```bash
python docker_output_to_timeline.py --docker_tsv_output_path <unsummarized Docker output TSV> --cancer_type <cancer type> --output_dir <output_dir>
```
Assuming successful completion the resulting summarized JSON can be passed to `eval_timeline.py` via the `--pred_path` parameter.

## Submission formats
### Checklists for formatting:
1. The file should include all patient ids in key of json file even if there is no predicted timeline.
2. Every key or values should be in lower case.
3. Saving in JSON format will automatically make all key/values in `str` type.
4. The evaluation script will need `--gold_path`, `--pred_path`, and `--all_id_path`. Note that `--gold_id_path` is reserved for test dataset - organizers will evaluate participants' submission with this option.
5. All timex values should be in [ISO 8601 standard format](https://en.wikipedia.org/wiki/ISO_8601), either YYYY-MM-DD and YYYY-Www(-dd), and **a week starts from Monday**. 
- In our gold annotated data, all values are normailzed with [clulab/timenorm](https://github.com/clulab/timenorm) library.

##### Note about ISO 8601 standard format (edge cases):
* In the relaxed_to_month setting, if a week spans two months, we consider both months as the answer span (e.g., 2024-W05 is from 2024-01-29 to 2024-02-04, so both January 2024 and February 2024 are correct).
* If there are unforeseen conflicts between the ISO 8601 standard and common sense regarding week representation (YYYY-Www(-dd)), we can manually review them and consider both answers to be correct. 


Example of article 1 (empty patient id case):
```python
{
    "patient01":[],
}
```
<hr>

## Installation

Clone this repository by
```bash
export EVAL_LIB_PATH=${HOME}
#export EVAL_LIB_PATH=<Or path to clone this repo>
cd ${EVAL_LIB_PATH}
git clone https://github.com/BCHHealthNLP/chemoTimelinesEval.git
``` 

### Environments
Please install required library by:
```bash
cd ${EVAL_LIB_PATH}/ChemoTimelines

pip install -r requirements.txt
```

We tested this evaluation codes on the following envirenments. 
* Ubuntu 22.04, python 3.10

<hr>

## Description

### Input:

- both gold and pred timelines are in json files, which is a dictionary maps patient name to the list of tuples of that patient, i.e. `{patient_name1: [tuple1, tuple2, ….], patient_name2, …}`
- tuple: `<chemotherapy_mention, temporal_relation, date>`, e.g. <carboplatin, contains-1, 2013-02-22>, we also call it <source, relation/label, target>.
- Note contains-1 is the inverse of contains, meaning carboplatin is contained by that date.

### Basic evaluation idea:

For each patient, we compare the tuples in gold timelines and in the predicted timelines, compute the F1 score.

### Setting

1. Strict evaluation
    1. a match means in a predicted tuple, the source, relation, and target are exactly the same as the ones in gold.
    2. For example, <carboplatin, contains-1, 2013-02-22> matches  <carboplatin, contains-1, 2013-02-22>, but doesn’t match  <carboplatin, **begins-on**, 2013-02-22> or  <carboplatin, contains-1, 2013-**03**-22> 
2. Relaxed evaluation:
    1. 3 settings: relaxed to date, relaxed to month, and relaxed to year
    2. for all three settings, we are relaxed about:
        1. **the label**, more specifically, we consider `contains-1` can be replaced by `begins-on` or `ends-on`, `begins-on` can be replaced by `contains-1`, `ends-on` can be replaced by `contains-1`, **but begins-on and ends-on cannot be replaced by each other.** 
        2. **the range**, meaning the predicted tuple fall in the correct range defined by the begins-on and ends-on dates in gold.
    3. for relaxed to month, we only care about if the predicted year-month matches the one in gold; for relaxed to year, we only care about if the predicted year matches the one in gold.
    
### Running the script
    
here is the example command to run evaluation on breast dev data.

```bash
export DATA_PATH=<path/to/data>
export PRED_PATH=<path/to/prediction>
export ID_PATH=<path/to/id files>
    
python eval_timeline.py \
    --gold_path ${DATA_PATH}/breast_dev_gold_timelines.json \
    --pred_path ${PRED_PATH}/breast_dev_system_timelines.json \
    --all_id_path ${ID_PATH}/breast_dev_all_ids.txt \
    --strict

# This option will show the official score
python eval_timeline.py \
    --gold_path ${DATA_PATH}/breast_dev_gold_timelines.json \
    --pred_path ${PRED_PATH}/breast_dev_system_timelines.json \
    --all_id_path ${ID_PATH}/breast_dev_all_ids.txt \
    --relaxed_to month
```
