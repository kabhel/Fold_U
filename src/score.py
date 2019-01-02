"""
    .. module:: Score
      :synopsis: This module implements the Score class.
"""

# Third-party modules
import os
import pandas as pd


def normalize_score(score_type):
    """
        Normalization of a score using the min-max scaling method (values between 0 and 1):

        .. math::

          zscore = {score - min(score) \over max(score) - min(score)}

        The Robust Z-Score for normalisation is not giving better results because it is usually
        used to better detect outliers, the formula is:

        .. math::

          { | x-med(x) | \over mad(x)}

        where:
        rz = Robust Zscore

        med(x) = Median

        mad = Median absolute deviation

        Args:
            score_type (Pandas Series): A score

        Returns:
            Pandas Series: The score normalized
    """
    try:
        normalized = (score_type - min(score_type)) / (max(score_type) - min(score_type))
    except ZeroDivisionError:
        print("\nError: normalize_score: Division by 0\n")
        normalized = 0
    return normalized
    # Use sklearn module instead
    # https://web.archive.org/web/20160520170701/http://chrisalbon.com:80/python/pandas_normalize_column.html


class Score:
    """
    .. class:: Score

      This class groups informations about a score.

    Attributes:
        iterator (iterator): An iterator of the generated scores.
    """

    def __init__(self, iterator):
        self.iterator = iterator

    def write_score(self, res_path, nb_pdb, alignment_dict):
        """
            Creates a scores.csv file containing for each line a template's name and its different
            scores. Creates the pdb files of the n best sum_scores.

            Args:
                res_path (str): The path of the directory where to stock the created files.
                nb_pdb (int): Number of pdb to create using the n first templates.
                alignment_dict (dictionary): A dictionary containing Alignment objects.
        """
        os.makedirs(res_path+"/pdb", exist_ok=True)

        # A dataframe is created with pandas and elements of the iterator are stored
        scores_df = pd.DataFrame(columns=['benchmark', 'alignment', 'threading', 'modeller',
                                          'secondary_structure', 'solvent_access', 'co_evolution'])
        for _, ali_score, thr_score, modeller_score, ss_score,\
                solvent_access_score, ccmpred_score, name, benchmark in sorted(self.iterator):
            scores_df.loc[name] = [benchmark, ali_score, thr_score, modeller_score,
                                   ss_score, solvent_access_score, ccmpred_score]

        # Normalization of the scores.
        # Not the ss_score neither solvent access_score because they are already between 0-1
        for index in ['alignment', 'threading', 'modeller', 'co_evolution']:
            scores_df[index] = normalize_score(scores_df[index])
        # Sum of the different scores and normalization
        scores_df['sum_scores'] = normalize_score(scores_df['alignment']
                                                  + scores_df['threading']
                                                  + scores_df['modeller']
                                                  + scores_df['secondary_structure']
                                                  + scores_df['solvent_access']
                                                  + scores_df['co_evolution'])

        # The machine learning was done using logistic regression with the R script located in
        # scripts/machine_learning.R
        # The script will find the best weights to apply to each score type in order to optimize
        # the benchmarking, that is to say it will give higher weights to scores which allow to
        # discriminate best the benchmarking proteins (of class Fold, Family and SuperFamily)
        # The Machine Learning script was run after the benchmarking results were generated.
        # The optimized weights are directly reported here to build a new column in CSV
        # scores.csv file for a weighted_combined_score.
        # The templates are then ranked according to this weighted_combined_scores.
        scores_df['weighted_combined_scores'] = normalize_score(
                                                    -10.8256
                                                  + 4.5026 * scores_df['alignment']
                                                  + -0.0764 * scores_df['threading']
                                                  + 1.0713 * scores_df['modeller']
                                                  + 3.8989 * scores_df['secondary_structure']
                                                  + -2.7788 * scores_df['solvent_access']
                                                  + -1.2846 * scores_df['co_evolution'])

        # Sort of the templates according to the weighted_combined_scores
        scores_df = scores_df.sort_values(by="weighted_combined_scores", ascending=False)
        # A csv file containing the normalized scores is created
        scores_df.to_csv(res_path+"/scores.csv")

        # Write the required ranking output with columns: Rank | Family | Score
        ranked_scores_df = pd.read_csv(res_path+"/scores.csv")
        ranks = pd.Series(list(range(1, len(ranked_scores_df)+1)))
        tmp = pd.concat([ranks, ranked_scores_df], axis=1)
        tmp.columns = ['Rank', 'Template', 'benchmark', 'alignment', 'threading', 'modeller',
                       'secondary_structure', 'solvent_access', 'co_evolution', 'sum_scores',
                       'weighted_combined_scores']
        ranked_scores_df = pd.concat([tmp["Rank"], tmp["Template"], tmp['weighted_combined_scores']], axis=1)
        with open(res_path+ "/ranking.txt", "w") as f_out:
            f_out.write("{}{:>22}{:>30}\n{}\n".format("Rank", "Template", "Weighted Combined Scores", "*"*56))
            for index, row in ranked_scores_df.iterrows():
                f_out.write("{:<3}{:>23}{:>12.4f}\n".format(str(row["Rank"]),
                                                            str(row["Template"]),
                                                            row["weighted_combined_scores"]))

        # Only nb_pdb pdb files are created
        for i in range(nb_pdb):
            pdb_filename = res_path + "/pdb/top_" + str(i+1) + ".pdb"
            alignment_dict[scores_df.index[i]].write_pdb(pdb_filename)
