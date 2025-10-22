"""
Hey. We're gonna load up the output of run_results_summary and make a nice LaTeX table from it.

"""

import os

import numpy as np
import pandas as pd

file_load_path = (
    "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes"
)
df = pd.read_hdf(os.path.join(file_load_path, "variability_results.h5"))

df.loc["n_nonperiodic_vars"] = df.loc["n_variables"] - df.loc["n_periodic_vars"]
df.loc["n_nonperiodic_vars_statistical"] = (
    df.loc["n_variables_statistical"] - df.loc["n_periodic_vars_statistical"]
)

df.loc["periodic_rate"] = (
    df.loc["n_periodic_vars_statistical"] / df.loc["n_statistical"]
)
df.loc["non_periodic_rate"] = (
    df.loc["n_nonperiodic_vars_statistical"] / df.loc["n_statistical"]
)
df.loc["variability_rate"] = df.loc["n_variables_statistical"] / df.loc["n_statistical"]


row_specs = [
    {"label": r"\textbf{Cluster membership}", "section": True},
    {"key": "n_members", "label": "Members identified by L16"},
    {
        "key": "n_M0_and_later",
        "label": r"Members with SpT $\geq$ M0",
        # "rule_after": True,
    },
    {
        "key": "n_IR_exc",
        "label": r"... and with mid-IR excess\tablenotemark{a}",
        "rule_after": True,
    },
    {"label": r"\textbf{Light curve quality and availability}", "section": True},
    {"key": "n_approved", "label": "Objects with acceptable light curves"},
    {
        "key": "n_statistical",
        "label": r"(subset with high-quality photometry)\tablenotemark{b}",
        "paren": True,
        "rule_after": True,
    },
    # {
    #     "key": "n_statistical",
    #     "label": r"(subset with high-quality photometry: $\sigma_K \leq 0.05$ mag, no error flags)",
    #     "paren": True,
    #     "rule_after": True,
    # },
    {"label": r"\textbf{Variability statistics}", "section": True},
    {
        "key": "n_automatic_vars",
        "label": r"Automatically determined variables",
        # "paren": True,
        # "rule_after": True,
    },
    {
        "key": "n_automatic_vars_statistical",
        "label": r"(high-quality)\tablenotemark{b}",
        "paren": True,
        "rule_after": True,
    },
    {
        "key": "n_subjective_vars",
        "label": r"Subjectively determined variables",
        # "paren": True,
        # "rule_after": True,
    },
    {
        "key": "n_subjective_vars_statistical",
        "label": r"(high-quality)\tablenotemark{b}",
        "paren": True,
        "rule_after": True,
    },
    {
        "key": "n_periodic_vars",
        "label": "Periodic variables",
        # "bold": True,
        # "rule_after": True,
    },
    {
        "key": "n_periodic_vars_statistical",
        "label": r"(high-quality)\tablenotemark{b}",
        "paren": True,
        # "bold": True,
        "rule_after": True,
    },
    # {
    #     "key": "n_nonperiodic_vars",
    #     "label": r"Non-periodic variables\tablenotemark{c}",
    #     "bold": True,
    #     # "rule_after": True,
    # },
    # {
    #     "key": "n_nonperiodic_vars_statistical",
    #     "label": r"(high-quality)\tablenotemark{b}",
    #     "paren": True,
    #     "bold": True,
    #     "rule_after": True,
    # },
    {
        "key": "n_variables",
        "label": r"All variables\tablenotemark{c}",
        # "paren": True,
        # "rule_after": True,
        "bold": True,
    },
    {
        "key": "n_variables_statistical",
        "label": r"(high-quality)\tablenotemark{b}",
        "paren": True,
        "rule_after": True,
        "bold": True,
    },
    {
        "label": r"\textbf{Variability rates (computed only on high-quality subset)}",
        "section": True,
    },
    {
        "key": "periodic_rate",
        "label": r"Periodic variability rate",
        # "paren": True,
        # "rule_after": True,
        # "bold": True,
        "rate": True,
    },
    {
        "key": "non_periodic_rate",
        "label": r"Non-periodic variability rate",
        # "paren": True,
        # "rule_after": True,
        # "bold": True,
        "rate": True,
    },    
    {
        "key": "variability_rate",
        "label": r"Total variability rate",
        # "paren": True,
        # "rule_after": True,
        "bold": True,
        "rate": True,
    },
]

end_matter = r"""
\tablenotetext{a}{
    Mid-infrared excess was assessed by L16 by inspecting \textit{Spitzer} data.
    An infrared excess is not used as a selection criterion for whether a source was
    included in our sample.
}
\tablenotetext{b}{
    We define an object to have ``high-quality photometry'' if both (a) its median $K$ 
    band uncertainty is less than 0.05 mag ($\sigma_K \leq 0.05$) and (b) all its data 
    are free of error processing flags.
}
\tablenotetext{c}{
    Many periodic variables were also identified as automatic variables; hence,
    the total variables count is smaller than the sum of the automatic, periodic, and 
    subjective variables counts.
}
\tablecomments{
    A full description of the variability identification procedure appears in \tsr{section}
}
"""


def to_deluxetable_custom(df, row_specs, title, label, end_matter):
    lines = []
    lines.append(r"\begin{deluxetable}{lccc}")
    lines.append(r"\tabletypesize{\scriptsize}")
    lines.append(r"\tablecaption{" + title + r"\label{" + label + r"}}")
    lines.append(
        r"\tablehead{ & \colhead{IC 348} & \colhead{NGC 1333} & \colhead{~~~Total~~~} }"
    )
    lines.append(r"\startdata")

    for spec in row_specs:
        # ---- Handle section headers ----
        if spec.get("section"):
            lines.append(r"\multicolumn{4}{l}{" + spec["label"] + r"} \\")
            lines.append(r"\hline")
            continue

        key = spec.get("key")
        # if key not in df.index:
        #     continue
        row = df.loc[key]

        # get ready
        ngc, ic, total = (
            (row["ngc"]),
            (row["ic"]),
            (row["total"]),
        )
        if spec.get('rate'):
            ngc, ic, total = [f"{v * 100:.1f}\\%" for v in (ngc, ic, total)]
        else:
            ngc, ic, total = [f"{int(v)}" for v in (ngc, ic, total)]

        # Parentheses for subset rows
        fmt = "({})" if spec.get("paren", False) else "{}"
        ngc, ic, total = (
            fmt.format(ngc),
            fmt.format(ic),
            fmt.format(total),
        )

        # Bold formatting
        if spec.get("bold", False):
            label_text = r"\textbf{" + spec["label"] + "}"
            ngc, ic, total = [r"\textbf{" + v + "}" for v in (ngc, ic, total)]
        else:
            label_text = spec["label"]

        # Data line
        lines.append(f"{label_text} & {ic} & {ngc} & {total} \\\\")

        # Optional horizontal rule
        if spec.get("rule_after", False):
            lines.append(r"\hline")

    lines.append(r"\enddata")
    lines.append(end_matter)
    lines.append(r"\end{deluxetable}")
    return "\n".join(lines)


if __name__ == "__main__":

    title = "Information on objects considered for this study"
    label = "tab:sample"

    lines = to_deluxetable_custom(df, row_specs, title, label, end_matter)
