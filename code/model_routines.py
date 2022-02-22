from tensorflow.keras import backend as K
from tensorflow.keras.metrics import SpecificityAtSensitivity
from tensorflow.keras.metrics import Precision
from tensorflow.keras.metrics import Recall, FalsePositives, FalseNegatives, TruePositives, TrueNegatives

def recall(y_true, y_pred):
    m = Recall()
    m.update_state(y_true, y_pred)
    pre_res = m.result()
    res = pre_res.numpy()
    return res


def f1(y_true, y_pred):
    m = Precision()
    m.update_state(y_true, y_pred)
    pre_precision = m.result()
    precision = pre_precision.numpy()

    m = Recall()
    m.update_state(y_true, y_pred)
    pre_recall = m.result()
    recall = pre_recall.numpy()

    return 2 * ((precision * recall) / (precision + recall + K.epsilon()))


def balanced_acc(y_true, y_pred):
    selectivity = recall(y_true, y_pred)
    specificity = SpecificityAtSensitivity(selectivity)
    specificity.update_state(y_true, y_pred)
    specificity = specificity.result().numpy()
    return (selectivity + specificity) / 2


def new_bac(y_true, y_pred):
    m = Recall()
    m.update_state(y_true, y_pred)
    pre_recall = m.result()
    recall = pre_recall.numpy()

    n = FalsePositives()
    n.update_state(y_true, y_pred)
    pre_res_2 = n.result()
    fp = pre_res_2.numpy()

    l = TrueNegatives()
    l.update_state(y_true, y_pred)
    pre_res = l.result()
    tn = pre_res.numpy()
    specificity = tn / (tn + fp)

    return (recall + specificity) / 2