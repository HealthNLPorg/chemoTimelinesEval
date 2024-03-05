import argparse
import json
import logging
import os
from collections import defaultdict
from datetime import datetime

import dateutil.parser

VERSION = "v20240223"

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
"""
    Set logger level as INFO for per patient-level evaluation results. 
    Set DEBUG for more detailed results.
"""


class Chemo:
    def __init__(self, text, first_start=None, last_end=None, cui=None):
        self.text = text
        self.first_start = first_start
        self.last_end = last_end
        self.cui = cui

    def __str__(self):
        return "\t".join(
            [self.text if self.text else "Null", self.cui if self.cui else "Null"]
        )


def read_all_patients(data_path):
    # Note that all key/value pairs of JSON are always of the type str.
    # https://docs.python.org/3/library/json.html
    with open(data_path, "r") as fr:
        all_patient_timelines = json.load(fr)
    return all_patient_timelines


def relaxed_rel_eval(incorrect, missing, preds, golds):
    not_truly_incorrect = []
    not_truly_missing = []
    for ptup in incorrect:
        is_not_truly_incorrect = False
        chemo, rel, timex = ptup
        # Basically we think contains-1 can be replaced by begins-on/ends-on,
        # and begins-on/ends-on can be replaced by contains-1.
        if rel in ["begins-on", "ends-on"]:
            if [chemo, "contains-1", timex] in golds:
                is_not_truly_incorrect = True
        elif rel == "contains-1":
            if [chemo, "begins-on", timex] in golds or [
                chemo,
                "ends-on",
                timex,
            ] in golds:
                is_not_truly_incorrect = True
        if is_not_truly_incorrect:
            not_truly_incorrect.append(ptup)

    for gtup in missing:
        is_not_truly_missing = False
        chemo, rel, timex = gtup

        if rel in ["begins-on", "ends-on"]:
            if [chemo, "contains-1", timex] in preds:
                is_not_truly_missing = True
        elif rel == "contains-1":
            if [chemo, "begins-on", timex] in preds or [
                chemo,
                "ends-on",
                timex,
            ] in preds:
                is_not_truly_missing = True
        if is_not_truly_missing:
            not_truly_missing.append(gtup)
    # print(len(not_truly_incorrect), len(not_truly_missing))
    return not_truly_incorrect, not_truly_missing


def relaxed_within_range_eval(incorrect, missing, gold_chemos, pred_chemos):
    """
        incorrect: false positive,
        missing: false negative,
        gold_chemos and pred_chemos, basically use Chemo object to get the start and end dates for each chemo
    """
    not_truly_incorrect = []
    not_truly_missing = []
    for ptup in incorrect:
        is_not_truly_incorrect = False
        source, rel, target = ptup

        target = target.replace("w", "W")
        if "W" in target:
            target = datetime.strptime(target + "-1", "%Y-W%W-%w")
        else:
            target = dateutil.parser.parse(target)

        if source in gold_chemos:
            gold_start, gold_end = (
                gold_chemos[source].first_start,
                gold_chemos[source].last_end,
            )
            if not gold_start or not gold_end:
                continue
            if rel in ["ENDS-ON", "ends-on"]:
                # The end date predicted by the system (target), is before the gold start date, then it's wrong.
                if target <= gold_start:
                    continue
            if rel in ["BEGINS-ON", "begins-on"]:
                # The start date predicted by the system (target) is after the gold end date, then it's wrong.
                if target >= gold_end:
                    continue
            # If the predicted date is in between the gold start and end date, i.e. in the correct range,
            # we consider this is correct, not truly false positive.
            if gold_start <= target <= gold_end:
                is_not_truly_incorrect = True
        if is_not_truly_incorrect:
            not_truly_incorrect.append(ptup)

    for gtup in missing:
        is_not_truly_missing = False
        source, rel, target = gtup

        target = target.replace("w", "W")
        if "W" in target:
            target = datetime.strptime(target + "-1", "%Y-W%W-%w")
        else:
            target = dateutil.parser.parse(target)

        if source in pred_chemos:
            pred_start, pred_end = (
                pred_chemos[source].first_start,
                pred_chemos[source].last_end,
            )
            if not pred_start or not pred_end:
                continue
            if rel in ["ENDS-ON", "ends-on"]:
                if target <= pred_start:
                    continue
            if rel in ["BEGINS-ON", "begins-on"]:
                if target >= pred_end:
                    continue
            # This is saying, for example, <taxol, contains-1, 2011-03-01> is missing in predictions, i.e. is false negative,
            # however, we can find <taxol, begins-on, 2011-01-01> and <taxol, ends-on, 2011-05-31> in predictions,
            # we consider <taxol, contains-1, 2011-03-01> is not false negative, because the system predicted the
            # correct range that covers the gold timeline.
            if pred_start <= target <= pred_end:
                is_not_truly_missing = True
        if is_not_truly_missing:
            not_truly_missing.append(gtup)
    return not_truly_incorrect, not_truly_missing


def group_chemo_dates(golds):
    gold_group_by_start_end = defaultdict(lambda: defaultdict(list))
    for tup in golds:
        source, label, target = tup
        if label.upper() not in ["BEGINS-ON", "ENDS-ON"]:
            continue
        target = target.replace("w", "W")
        if "-W" in target:
            target = datetime.strptime(target + "-1", "%Y-W%W-%w")
        else:
            target = dateutil.parser.parse(target)

        gold_group_by_start_end[source][label].append(target)
    all_gold_chemos = {}
    # all_gold_chemos: maps from the text of this chemo to the chemo event obj
    for chemo, labels in gold_group_by_start_end.items():
        first_date = (
            min(labels["BEGINS-ON".lower()]) if "BEGINS-ON".lower() in labels else None
        )
        last_date = (
            max(labels["ENDS-ON".lower()]) if "ENDS-ON".lower() in labels else None
        )
        # For each chemo, find the earliest start date, and the latest end date,
        # use them to get the span of this chemo, the span will be used when doing relaxed evaluation.
        chemo_event = Chemo(text=chemo, first_start=first_date, last_end=last_date)
        all_gold_chemos[chemo] = chemo_event
    return all_gold_chemos


def normalize_to_month_and_year(golds):
    month_only_pairs = []
    year_only_pairs = []
    for tup in golds:
        source, label, target = tup
        target = target.replace("w", "W")
        if "W" in target:
            target = datetime.strptime(target + "-1", "%Y-W%W-%w")
        else:
            target = dateutil.parser.parse(target)
        year = target.year
        month = target.month
        if month < 10:
            normalized_month = str(year) + "-0" + str(month)
        else:
            normalized_month = str(year) + "-" + str(month)
        normalized_year = str(year)

        month_pair = [source, label, normalized_month]
        year_pair = [source, label, normalized_year]

        if month_pair not in month_only_pairs:
            month_only_pairs.append(month_pair)
        if year_pair not in year_only_pairs:
            year_only_pairs.append(year_pair)

    has_more_specific_chemos_month = summarize_timeline(month_only_pairs)
    has_more_specific_chemos_year = summarize_timeline(year_only_pairs)
    month_only_pairs = [
        tup for tup in month_only_pairs if tup not in has_more_specific_chemos_month
    ]
    year_only_pairs = [
        tup for tup in year_only_pairs if tup not in has_more_specific_chemos_year
    ]
    return month_only_pairs, year_only_pairs


def summarize_timeline(timelines):
    """
        This is to postprocess timelines one more time after we normalized original timeline
        to year only or month only timelines. What it does is: if we have a generic chemo mention,
        e.g. <chemotherapy, contains-1, 2011-01>, or <chemoradiation, contains-1, 2011-01>, we want to see
        if we can have more specific chemo mention happened on the same date with the same label,
        e.g. <Taxol, contains-1, 2011-01>. If we find a more specific chemo mention,
        we would ignore the generic chemo mention, only add <Taxol, contains-1, 2011-01> to the timeline.
    """
    date_rel_to_chemo = defaultdict(lambda: defaultdict(list))
    for tup in timelines:
        chemo, rel, date = tup
        date_rel_to_chemo[date][rel].append(chemo)

    has_more_specific_chemos = []
    for date, rel_chemos in date_rel_to_chemo.items():
        for rel, chemos in rel_chemos.items():
            for chemo in chemos:
                # chemo.startswith("chemo") is how we check if this is a generic chemo mention
                if chemo.startswith("chemo"):
                    if len(date_rel_to_chemo[date][rel]) > 1:
                        has_more_specific_chemos.append([chemo, rel, date])
    return has_more_specific_chemos


def compute_f(true_pos, false_pos, false_neg):
    precision = len(true_pos) / (len(true_pos) + len(false_pos))
    recall = len(true_pos) / (len(true_pos) + len(false_neg))
    f1 = 2 * (precision * recall) / (precision + recall)
    return f1


def strict_eval(gold, pred):
    true_pos = [prediction for prediction in pred if prediction in gold]
    false_pos = [prediction for prediction in pred if prediction not in gold]
    false_neg = [correct for correct in gold if correct not in pred]
    not_truly_fp = fp_fn_single_count(false_pos, false_neg)
    false_pos = [pred for pred in false_pos if pred not in not_truly_fp]
    return true_pos, false_pos, false_neg


def fp_fn_single_count(false_pos, false_neg):
    """
        What it does here is: let's say in pred we have <Taxol, BEGINS-ON, 2011-01-01>,
        in gold we have <Taxol, CONTAINS-1, 2011-01-01>, then <Taxol, BEGINS-ON, 2011-01-01> would be false positive,
        <Taxol, CONTAINS-1, 2011-01-01> would be false negative, that means, the same mistake is counted twice,
        once in fp, once in fn. So, here, we want to make sure, we count <Taxol, CONTAINS-1, 2011-01-01> as false negative,
        and don't count <Taxol, BEGINS-ON, 2011-01-01> as false positive.
    """
    not_truly_fp = []
    # false_neg_tracker: (chemo, timex) to label
    false_neg_tracker = {(item[0], item[-1]): item[1] for item in false_neg}
    for ptup in false_pos:
        # E.g. we check in <Taxol, 2011-01-01> is already in false negative.
        if (ptup[0], ptup[-1]) in false_neg_tracker:
            not_truly_fp.append(ptup)
    return not_truly_fp


def relaxed_eval(gold, gold_chemo, pred, pred_chemo):
    true_pos = [prediction for prediction in pred if prediction in gold]
    false_pos = [prediction for prediction in pred if prediction not in gold]
    false_neg = [correct for correct in gold if correct not in pred]
    not_truly_fp_with_range, not_truly_fn_with_range = relaxed_within_range_eval(
        incorrect=false_pos,
        missing=false_neg,
        gold_chemos=gold_chemo,
        pred_chemos=pred_chemo,
    )
    not_truly_fp_with_label, not_truly_fn_with_label = relaxed_rel_eval(
        incorrect=false_pos, missing=false_neg, preds=pred, golds=gold
    )

    not_truly_fp_as_label_single_count = fp_fn_single_count(false_pos, false_neg)

    truly_tp = true_pos
    truly_fp, truly_fn = [], []
    for tup in false_pos:
        if tup in not_truly_fp_with_range or tup in not_truly_fp_with_label:
            # Add this one to true positive if it's not considered as true fp
            truly_tp.append(tup)
        elif tup in not_truly_fp_as_label_single_count:
            continue
        else:
            truly_fp.append(tup)
    for tup in false_neg:
        if tup in not_truly_fn_with_range or tup in not_truly_fn_with_label:
            continue
        truly_fn.append(tup)
    return (
        truly_tp,
        truly_fp,
        truly_fn,
        not_truly_fp_with_range,
        not_truly_fn_with_range,
        not_truly_fp_with_label,
        not_truly_fn_with_label,
    )


def evaluation(gold, pred, args):
    # Get the earliest start and latest end dates for each chemo
    all_gold_chemos = group_chemo_dates(gold)
    all_pred_chemos = group_chemo_dates(pred)

    gold_month_timeline, gold_year_timeline = normalize_to_month_and_year(gold)
    pred_month_timeline, pred_year_timeline = normalize_to_month_and_year(pred)

    if args.strict:
        true_pos, false_pos, false_neg = strict_eval(gold, pred)
    else:
        assert (
            args.relaxed_to
        ), "For relaxed evaluation, please specify --relaxed_to flag"
        if args.relaxed_to == "day":
            # rmv_from_fp_range: remove from false positive because it's in the right range;
            # rmv_from_fp_label: remove from false positive because the label is consider correct.
            (
                true_pos,
                false_pos,
                false_neg,
                rmv_from_fp_range,
                rmv_from_fn_range,
                rmv_from_fp_label,
                rmv_from_fn_label,
            ) = relaxed_eval(gold, all_gold_chemos, pred, all_pred_chemos)
        elif args.relaxed_to == "month":
            (
                true_pos,
                false_pos,
                false_neg,
                rmv_from_fp_range,
                rmv_from_fn_range,
                rmv_from_fp_label,
                rmv_from_fn_label,
            ) = relaxed_eval(
                gold_month_timeline,
                all_gold_chemos,
                pred_month_timeline,
                all_pred_chemos,
            )
        elif args.relaxed_to == "year":
            (
                true_pos,
                false_pos,
                false_neg,
                rmv_from_fp_range,
                rmv_from_fn_range,
                rmv_from_fp_label,
                rmv_from_fn_label,
            ) = relaxed_eval(
                gold_year_timeline, all_gold_chemos, pred_year_timeline, all_pred_chemos
            )
        else:
            raise ValueError("--relaxed_to must be one of 'day', 'month', or 'year'")

    logger.debug(f"true_pos... {len(true_pos)}")
    for item in true_pos:
        logger.debug(f"{item}")
        pass
    logger.debug(f"false_pos... {len(false_pos)}")
    for item in false_pos:
        logger.debug(f"{item}")
        pass
    logger.debug(f"false_neg... {len(false_neg)}")
    for item in false_neg:
        logger.debug(f"{item}")
        pass

    if not args.strict:
        logger.debug(f"removed from false_pos_range... {len(rmv_from_fp_range)}")
        logger.debug(f"{rmv_from_fp_range}")
        logger.debug(f"removed from false_pos_label... {len(rmv_from_fp_label)}")
        logger.debug(f"{rmv_from_fp_label}")

        logger.debug(f"removed from false_neg_range... {len(rmv_from_fn_range)}")
        logger.debug(f"{rmv_from_fn_range}")
        logger.debug(f"removed from false_neg_label... {len(rmv_from_fn_label)}")
        logger.debug(f"{rmv_from_fn_label}")

    if len(true_pos) + len(false_neg) == 0:
        precision, recall, f1 = 0, 0, 0
    elif len(true_pos) + len(false_pos) == 0:
        precision, recall, f1 = 0, 0, 0
    else:
        precision = len(true_pos) / (len(true_pos) + len(false_pos))
        recall = len(true_pos) / (len(true_pos) + len(false_neg))
        if precision + recall:
            f1 = 2 * (precision * recall) / (precision + recall)
        else:
            f1 = 0

    logger.info(f"precision: {precision}")
    logger.info(f"recall: {recall}")
    logger.info(f"f1: {f1}")

    return true_pos, false_pos, false_neg, precision, recall, f1


def read_files(args):
    gold_all_patient = read_all_patients(args.gold_path)
    pred_all_patient = read_all_patients(args.pred_path)

    with open(args.all_id_path) as fp:
        all_ids = [line.splitlines()[0] for line in fp.readlines()]

    # Sanity check
    if len(all_ids) == 0 or len(all_ids) != len(pred_all_patient):
        raise ValueError(
            "Malformated or some patients are missing in prediction file."
        )
    if not args.gold_id_path:
        if len(all_ids) != len(gold_all_patient):
            raise ValueError(
                "Error in gold annotated ids file. Check the content of file in gold_id_path. Length of all_ids: %s, gold_all_patient: %s"
                % (len(all_ids), len(gold_all_patient))
            )

    # Only for orgnizers: for test dataset - screening silver datasets
    if args.gold_id_path:
        if not os.path.exists(args.gold_id_path):
            raise ValueError("Error in gold annotated ids file path")
        with open(args.gold_id_path) as fp:
            gold_ids = [line.splitlines()[0] for line in fp.readlines()]

        if (
            len(gold_ids) == 0
            or len(gold_ids) != len(gold_all_patient)
            or len(gold_ids) > len(pred_all_patient)
        ):
            raise ValueError("Error in gold annotated ids file. Check the content of file in gold_id_path")
    else:
        gold_ids = list(gold_all_patient.keys())

    return pred_all_patient, gold_all_patient, all_ids, gold_ids


def micro_average_metrics(
        all_true_pos:dict, 
        all_false_pos:dict, 
        all_false_neg:dict
        ) -> float:
    # Micro average metrics
    logger.info(
        f"tp, fp, fn over all patients: {sum(all_true_pos.values())}, {sum(all_false_pos.values())}, {sum(all_false_neg.values())}"
    )
    
    if len(all_true_pos) + len(all_false_pos) != 0:
        micro_precision = sum(all_true_pos.values()) / (sum(all_true_pos.values()) + sum(all_false_pos.values()))
    else:
        micro_precision = 0
    if len(all_true_pos) + len(all_false_neg) != 0:
        micro_recall = sum(all_true_pos.values()) / (sum(all_true_pos.values()) + sum(all_false_neg.values()))
    else:
        micro_recall = 0
    if micro_precision + micro_recall:
        micro_f1 = 2 * (micro_precision * micro_recall) / (micro_precision + micro_recall)
    else:
        micro_f1 = 0

    print("Micro average metrics")
    print("Micro precision:", micro_precision)
    print("Micro recall:", micro_recall)
    print("Micro f1:", micro_f1)
    print()

    return micro_f1


def macro_average_metrics(
        local_precision:dict, 
        local_recall:dict, 
        local_f1:dict,
        ) -> tuple:
    # Macro average metrics
    print("Macro average metrics")

    print("Type A evaluation: including the notes with no true relations")
    total_patients = len(local_f1)
    assert len(local_precision) == len(local_recall) == len(local_f1) == total_patients, \
        "total_patients should be equal to the number of local metrics"
    
    type_a_macro_prec = sum(local_precision.values()) / total_patients
    type_a_macro_rec = sum(local_recall.values()) / total_patients
    type_a_macro_f1 = sum(local_f1.values()) / total_patients
    print("[Type A] Macro precision:", type_a_macro_prec)
    print("[Type A] Macro recall:", type_a_macro_rec)
    print("[Type A] Macro f1: ", type_a_macro_f1)
    print()

    print("Type B evaluation: excluding the notes with no true relations")
    type_b_prec = [
        score for pat_id, score in local_precision.items() if local_relations[pat_id] != 0
    ]
    type_b_rec = [
        score for pat_id, score in local_recall.items() if local_relations[pat_id] != 0
    ]
    type_b_f1 = [
        score for pat_id, score in local_f1.items() if local_relations[pat_id] != 0
    ]
    assert len(type_b_prec) == len(type_b_rec) == len(type_b_f1), \
        "The number of local metrics should be the same."
    
    type_b_macro_prec = sum(type_b_prec) / len(type_b_prec)
    type_b_macro_rec = sum(type_b_rec) / len(type_b_rec)
    type_b_macro_f1 = sum(type_b_f1) / len(type_b_f1)
    print("[Type B] Macro precision:", type_b_macro_prec)
    print("[Type B] Macro recall:", type_b_macro_rec)
    print("[Type B] Macro f1: ", type_b_macro_f1)
    print()

    return type_a_macro_f1, type_b_macro_f1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Evaluate predicted output against gold annotations"
    )

    parser.add_argument(
        "--gold_path", required=True, help="A gold annotation json file"
    )
    parser.add_argument(
        "--pred_path", required=True, help="A predicted output json file"
    )
    parser.add_argument(
        "--all_id_path",
        required=True,
        help="Path to file with list of ids, delimited by new line characters",
    )
    parser.add_argument(
        "--gold_id_path",
        required=False,
        help="(Only for test evaluation) Path to file with list of gold annotated ids, delimited by new line characters",
    )

    parser.add_argument("--strict", action="store_true", help="do strict eval")
    parser.add_argument(
        "--relaxed_to",
        help="Type 'year' to only evaluate year, 'month' to evaluate year and month, "
        "or 'day' to evaluate year-month-day",
        choices=["day", "month", "year"],
    )
    args = parser.parse_args()

    logger.info(args)

    print(f"Evaluation code for ChemoTimelines Shared Task. Version: {VERSION}")
    print(f"Reading from files...")
    pred_all_patient, gold_all_patient, all_ids, gold_ids = read_files(args=args)
    print(f"predicted output: {len(pred_all_patient)}, gold annotation: {len(gold_all_patient)}, all ids: {len(all_ids)}")
    print()

    all_true_pos, all_false_pos, all_false_neg = {}, {}, {}
    local_relations = {} # Key = patient ID, Value = number of timeline in the patient
    local_f1, local_precision, local_recall = {}, {}, {} # Key = patient ID, Value = score for the patient

    for pred_patient, pred_timeline in pred_all_patient.items():
        # pred_patient: patient ID; pred_timeline: a list of <chemo, label, timex>
        if pred_patient not in gold_ids:
            continue

        if pred_patient not in gold_all_patient:
            raise ValueError(
                "The given patient ID: '%s' does not exist in the gold annotated dataset."
                % pred_patient
            )
        gold_timeline = gold_all_patient.get(pred_patient, [])

        # Handling patients without timeline. Models are expected to output "no timeline".
        """
            For type A evaluation, metric includes the patients with no gold timelines:
            If models successfully make no prediction, they will achieve scores of 
            1.0, 1.0, and 1.0 for the given patient in local precision, recall, and F1, respectively.
        """
        if len(gold_timeline) == 0:
            if len(pred_timeline) == 0:
                true_pos, false_pos, false_neg = [], [], []
                p, r, f_score = 1, 1, 1
            else:
                true_pos = []
                false_pos = pred_timeline
                false_neg = []
                p, r, f_score = 0, 0, 0
            local_relations[pred_patient] = 0
        else:
            logger.info(f"pred_patient ID: {pred_patient}")
            true_pos, false_pos, false_neg, p, r, f_score = evaluation(
                gold=gold_timeline, pred=pred_timeline, args=args
            )
            local_relations[pred_patient] = len(gold_timeline)

        all_true_pos[pred_patient] = len(true_pos)
        all_false_pos[pred_patient] = len(false_pos)
        all_false_neg[pred_patient] = len(false_neg)

        local_precision[pred_patient] = p
        local_recall[pred_patient] = r
        local_f1[pred_patient] = f_score


    _ = micro_average_metrics(
        all_true_pos=all_true_pos, 
        all_false_pos=all_false_pos, 
        all_false_neg=all_false_neg
    )

    type_a_macro_f1, type_b_macro_f1 = macro_average_metrics(
        local_precision=local_precision, 
        local_recall=local_recall, 
        local_f1=local_f1,
    )

    print("Official Score: Arithmetic mean of two types of Macro F1, type A and B, " + \
          "in 'relaxed to month' setting will be used for the leaderboard. ")
    if not (args.strict) and args.relaxed_to == "month":
        print("Official Score: ", (type_a_macro_f1 + type_b_macro_f1) / 2)
    else:
        print("To see the official score, please run without --strict flag, and set 'relaxed to month' setting by --relaxed_to=month ")

    print("Evaluation completed!")
