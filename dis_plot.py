

import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.special import boxcox1p
import pickle


ref = 'AAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCCGCCTGGCATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGGTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAGCGGTTTGACTCACGGGGATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCAAAATCAACGGGACTTTCCAAAATGTCGTAACAACTCCGCCCCATTGACGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCTCTCTGGCTAACTAGAGAACCCACTGCTTACTGGCTTATCGAAATTAATACGACTCACTATAGGGAGACCCAAGCTGGCTAGCGTTTAAACTTAAGCTTGGTACCGAGCTCGGATCCACTAGTCCAGTGTGGTGGAATTCTGCAGATATCCAGCACAGTGGCGGCCGCGCGCGGGCCCAGCCGGCGCCCGTCCGCGCCGGGCCCGTTTAAACCCGCTGATCAGCCTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTCTGAGTAGGTGTCATTCTATTCTGGGGGGTGGGGTGGGGCAGGACAGCAAGGGGGAGGATTGGGAAGACAATAGCAGGCATGCTGGGGATGCGGTGGGCTCTATGGCTTCTGAGGCGGAAAGAACCAGCTGGGGCTCTAGGGGGTATCCCCACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACCCCAAAAAACTTGATTAGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTAATTCTGTGGAATGTGTGTCAGTTAGGGTGTGGAAAGTCCCCAGGCTCCCCAGCAGGCAGAAGTATGCAAAGCATGCATCTCAATTAGTCAGCAACCAGGTGTGGAAA'

def calculate_signal_diff(adduct_signal, control_signal, pos_list, p_thresh):
    """
    计算每个碱基位点的：平均信号、倍数变化（FC）、t检验p值、显著性标记
    :param adduct_signal: 加合物去噪信号矩阵
    :param control_signal: 对照去噪信号矩阵
    :param pos_list: 目标位点列表
    :param p_thresh: p值阈值
    :return: 差异统计结果DataFrame
    """
    diff_result = []
    for i, pos in enumerate(pos_list):
        # 提取单个位点的所有reads信号
        a_sig = adduct_signal[i]
        c_sig = control_signal[i]
        # 计算均值
        a_mean = np.mean(a_sig)
        c_mean = np.mean(c_sig)
        # 倍数变化（FC=加合物均值/对照均值）
        fc = a_mean / c_mean if c_mean != 0 else np.inf
        # 两独立样本t检验（检验信号差异是否显著）
        t_stat, p_val = stats.ttest_ind(a_sig, c_sig, equal_var=False)  # 非等方差t检验
        # 显著性标记
        significant = "Yes" if p_val < p_thresh else "No"

        # 计算标准化偏差
        c_std = np.std(c_sig)
        mean_diff = np.mean(np.abs(a_sig - c_mean)/c_std)
        c_mean_diff = np.mean(np.abs(c_sig - c_mean)/c_std)
        diff_result.append({
            "position": pos,
            "adduct_mean": a_mean,
            "control_mean": c_mean,
            "fold_change": fc,
            "t_stat": t_stat,
            "p_value": p_val,
            "significant": significant,
            "mean_diff": mean_diff,
            "c_mean_diff": c_mean_diff,
        })
    # 转换为DataFrame并保存
    diff_df = pd.DataFrame(diff_result)
    return diff_df

# 执行信号差异统计
def calculate_pi_value(diff_df):
    """
    计算π值（DIS核心）：整合倍数变化（FC）和t检验p值，量化信号差异的强度+可靠性
    公式参考Xiao et al. (2014)：π = log2(FC) * (-log10(p_val))，符号反映差异方向
    :param diff_df: 信号差异统计DataFrame
    :return: 含π值的DIS结果DataFrame
    """
    dis_df = diff_df.copy()
    # 处理极端值（FC=0/inf，p_val=0）
    dis_df["log2_fc"] = np.log2(dis_df["fold_change"].replace([0, np.inf], [1e-6, 1e6]))
    dis_df["-log10_p"] = -np.log10(dis_df["p_value"].replace(0, 1e-300))
    # 计算π值（DIS核心指标，即每个碱基位点的DIS值）
    # dis_df["pi_value"] = dis_df["log2_fc"] * dis_df["-log10_p"]
    dis_df["pi_value"] = dis_df["-log10_p"]
    dis_df["mean_dff"] = dis_df["mean_diff"]
    dis_df["c_mean_dff"] = dis_df["c_mean_diff"]
    return dis_df


def plot_dis_radar2(dis_df, output_dir):
    """
    :param dis_df: 含pi_value的DIS结果DataFrame
    :param output_dir: 图片输出路径
    """
    # 准备数据：按位点排序，提取pi值
    dis_df = dis_df.sort_values("position").reset_index(drop=True)
    positions = [str(p) for p in dis_df["position"]]
    pi_values = dis_df["mean_dff"].values
    # 雷达图角度设置
    angles = np.linspace(-1/4 * np.pi, 5/4 * np.pi,len(positions), endpoint=True)[::-1]
    angles = np.where(angles < 0, angles + 2 * np.pi, angles)

    _, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    # 绘制DIS曲线
    ax.plot(angles, pi_values, "o-", linewidth=2, color="#e74c3c")
    ax.set_ylim(0, max(pi_values) * 1.2)
    ax.plot(angles[positions.index('506')], pi_values[positions.index('506')], '*', markersize=15, color='red')

    # 绘制0值参考线
    pi_values = dis_df["c_mean_dff"].values
    ax.plot(angles, pi_values, "o-", linewidth=2, color="k")

    labels = [f"{int(p)}:{ref[int(p)-1]}" for p in positions]
    ax.set_xticks(angles)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_rlabel_position(angles[0]*180/np.pi+5)
    ax.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, "DIS_radar_plot.png"), 
                dpi=300, bbox_inches="tight")
    # plt.show()
    plt.close()
    print(f"DIS雷达图保存至：{os.path.join(output_dir, 'DIS_radar_plot.png')}")


def denoise_signal(signal_mat):
    boxcox_lmbda = stats.boxcox_normmax(signal_mat + 1, method="mle")
    denoised_signal = boxcox1p(signal_mat + 1, boxcox_lmbda)
    return denoised_signal


if __name__ == "__main__":
    
    platform="QT"
    if platform=='ONT':
        start_position=11
    else:
        start_position=7
    root_dir = os.path.dirname(os.path.abspath(__file__))
    fig_save_dir = os.path.join(root_dir, 'fig', platform)
    os.makedirs(fig_save_dir, exist_ok=True)

    control_res_path = os.path.join(root_dir, 'data', f'{platform}_control_align_res.pkl')
    dAAL_res_path = os.path.join(root_dir, 'data', f'{platform}_dAAL_align_res.pkl')

    with open(control_res_path, 'rb') as f:
        control_res = pickle.load(f)['sigs']
        control_res = np.array(control_res)
    with open(dAAL_res_path, 'rb') as f:
        dAAL_res = pickle.load(f)['sigs']
        dAAL_res = np.array(dAAL_res)

    min_pos, max_pos = 475, 515
    resample = 6
    all_dAAL_sigs = []
    all_control_sigs = []

    for j in range(min_pos, max_pos):
        
        # 提取信号
        control_sigs = control_res[:, j*resample:(j+1)*resample].flatten()
        dAAL_sigs = dAAL_res[:, j*resample:(j+1)*resample].flatten()

        all_dAAL_sigs.append(dAAL_sigs)
        all_control_sigs.append(control_sigs)

    diff_df = calculate_signal_diff(all_dAAL_sigs, all_control_sigs, list(np.arange(min_pos, max_pos).astype(int)+12), p_thresh=0.05)
    dis_df = calculate_pi_value(diff_df)
    plot_dis_radar2(dis_df, fig_save_dir)