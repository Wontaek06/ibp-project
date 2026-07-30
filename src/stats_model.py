"""
Group comparison stats + a small classifier.

Sample sizes here are small (n=13-18 species per full dataset, fewer per
AFP-type group), so:
  - Kruskal-Wallis (non-parametric) is used instead of ANOVA
  - Dunn's post-hoc test needs `scikit-posthocs` (pip install scikit-posthocs)
  - the RandomForest classifier is evaluated with Leave-One-Out CV and is
    meant to surface FEATURE IMPORTANCE (what best explains AFP type),
    not to be reported as a high-accuracy predictive model - state this
    explicitly in the report given the sample size.
"""
import numpy as np
from scipy.stats import kruskal
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.metrics import accuracy_score


def group_comparison(df, value_col, group_col="afp_type"):
    """Kruskal-Wallis across groups for one variable. Returns None if <2 groups have data."""
    d = df.dropna(subset=[value_col, group_col])
    groups = [d[d[group_col] == t][value_col].values for t in d[group_col].unique()]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2:
        return None
    stat, p = kruskal(*groups)
    return {"variable": value_col, "H": round(float(stat), 2), "p": round(float(p), 4),
            "n_groups": len(groups), "n_total": len(d)}


def bootstrap_ci(values, statistic=np.median, n_boot=2000, ci=95, seed=42):
    """
    Percentile bootstrap CI for a group statistic (median by default).
    With n as small as 3-4 per group the interval will be wide - that is
    the honest result and should be shown rather than hidden.
    Returns (low, high), or (None, None) if fewer than 3 values.
    """
    vals = np.asarray([v for v in values if v is not None and not np.isnan(v)], dtype=float)
    if len(vals) < 3:
        return None, None
    rng = np.random.default_rng(seed)
    boots = [statistic(rng.choice(vals, size=len(vals), replace=True)) for _ in range(n_boot)]
    lo = np.percentile(boots, (100 - ci) / 2)
    hi = np.percentile(boots, 100 - (100 - ci) / 2)
    return round(float(lo), 3), round(float(hi), 3)


def posthoc_dunn(df, value_col, group_col="afp_type"):
    """Pairwise Dunn's test (Holm-corrected). Returns None if scikit-posthocs unavailable."""
    try:
        import scikit_posthocs as sp
    except ImportError:
        print("scikit-posthocs not installed - run: pip install scikit-posthocs")
        return None
    d = df.dropna(subset=[value_col, group_col])
    if d[group_col].nunique() < 2:
        return None
    return sp.posthoc_dunn(d, val_col=value_col, group_col=group_col, p_adjust="holm")


def classify_afp_type(df, feature_cols, target_col="afp_type"):
    """
    RandomForest + Leave-One-Out CV over feature_cols -> target_col.
    Returns (loo_accuracy, feature_importances_dict) or (None, None) if
    too few complete rows to fit meaningfully (< 6).
    """
    d = df.dropna(subset=feature_cols + [target_col])
    if len(d) < 6:
        return None, None

    X = d[feature_cols].values
    y = d[target_col].values

    clf = RandomForestClassifier(n_estimators=300, random_state=42)
    y_pred = cross_val_predict(clf, X, y, cv=LeaveOneOut())
    acc = accuracy_score(y, y_pred)

    clf.fit(X, y)  # refit on all data for feature importance
    importances = dict(zip(feature_cols, np.round(clf.feature_importances_, 3)))

    return round(float(acc), 3), importances
