

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle


ref = 'AAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCCGCCTGGCATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGGTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAGCGGTTTGACTCACGGGGATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCAAAATCAACGGGACTTTCCAAAATGTCGTAACAACTCCGCCCCATTGACGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCTCTCTGGCTAACTAGAGAACCCACTGCTTACTGGCTTATCGAAATTAATACGACTCACTATAGGGAGACCCAAGCTGGCTAGCGTTTAAACTTAAGCTTGGTACCGAGCTCGGATCCACTAGTCCAGTGTGGTGGAATTCTGCAGATATCCAGCACAGTGGCGGCCGCGCGCGGGCCCAGCCGGCGCCCGTCCGCGCCGGGCCCGTTTAAACCCGCTGATCAGCCTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTCTGAGTAGGTGTCATTCTATTCTGGGGGGTGGGGTGGGGCAGGACAGCAAGGGGGAGGATTGGGAAGACAATAGCAGGCATGCTGGGGATGCGGTGGGCTCTATGGCTTCTGAGGCGGAAAGAACCAGCTGGGGCTCTAGGGGGTATCCCCACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACCCCAAAAAACTTGATTAGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTAATTCTGTGGAATGTGTGTCAGTTAGGGTGTGGAAAGTCCCCAGGCTCCCCAGCAGGCAGAAGTATGCAAAGCATGCATCTCAATTAGTCAGCAACCAGGTGTGGAAA'


if __name__ == "__main__":
    z_thres_mix = 2.0

    for platform in ['QT', 'ONT']:
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


        resample = 6
        min_pos = 20
        max_pos = 1320

        c_means = np.zeros(max_pos-min_pos)
        c_stds = np.ones(max_pos-min_pos)

        control_zscores = np.zeros((len(control_res), max_pos-min_pos))
        dAAL_zscores = np.zeros((len(dAAL_res), max_pos-min_pos))
        for j in range(min_pos, max_pos):
            
            # 提取信号
            control_sigs = control_res[:len(control_res)//2, j*resample:(j+1)*resample].flatten()
            c_means[j-min_pos] = np.mean(control_sigs)
            c_stds[j-min_pos] = np.std(control_sigs)

            for i in range(len(control_res)):
                control_zscores[i, j-min_pos] = np.abs(np.mean(control_res[i, j*resample:(j+1)*resample]) - c_means[j-min_pos]) / c_stds[j-min_pos]

            for i in range(len(dAAL_res)):
                dAAL_zscores[i, j-min_pos] = np.abs(np.mean(dAAL_res[i, j*resample:(j+1)*resample]) - c_means[j-min_pos]) / c_stds[j-min_pos]

        threshs = [1, 2, 3]

        for z_thres in threshs:
            plt.figure(figsize=(10, 8))
            plt.title('zscore threshold: %.1f' % z_thres)
            out_ratio = np.zeros(max_pos-min_pos)
            for j in range(min_pos, max_pos):
                out_ratio[j-min_pos] = np.sum(control_zscores[len(control_zscores)//2:, j-min_pos] > z_thres) / (len(control_zscores)//2)
            
            plt.plot(np.arange(min_pos, max_pos)+start_position+1, out_ratio, color='k', label='control')

            for j in range(min_pos, max_pos):
                out_ratio[j-min_pos] = np.sum(dAAL_zscores[:, j-min_pos] > z_thres) / dAAL_zscores.shape[0]
            # # out_ratio = np.convolve(out_ratio, np.ones(5)/5, mode='same')
            plt.plot(np.arange(min_pos, max_pos)+start_position+1, out_ratio, color='r', label='dA-AL-II')
            plt.legend()
            plt.xlabel('Position')
            plt.ylabel('Ratio')
            plt.savefig(os.path.join(fig_save_dir, 'zscore_threshold_%d_out_ratio.png' % z_thres), dpi=300)
            plt.close()
        control_feas = control_zscores[len(control_zscores)//2:, 506-1-start_position-min_pos]
        dAAL_feas = dAAL_zscores[:, 506-1-start_position-min_pos]

        # 防止两组数据数量相差过大
        if len(control_feas) > len(dAAL_feas) * 1.2:
            control_feas = control_feas[:int(len(dAAL_feas)*1.2)]

        feas = np.concatenate([control_feas, dAAL_feas], axis=0)
        y_true = np.concatenate([np.zeros(len(control_feas)), np.ones(len(dAAL_feas))], axis=0)

        from sklearn.metrics import (confusion_matrix, ConfusionMatrixDisplay, auc, roc_curve)
        for z_thres in threshs:
            # plt.figure(figsize=(10, 8))
            pred_labels = np.array(feas >= z_thres).astype(bool)
            cm = confusion_matrix(y_true, pred_labels)
            tn, fp, fn, tp = cm.ravel()
            # 核心指标计算
            fpr = fp / (fp + tn)          # 假阳性率 (False Positive Rate)
            tpr = tp / (tp + fn)          # 真阳性率/召回率 (Recall)
            precision = tp / (tp + fp)    # 精确率 (Precision)
            f1_score = 2 * (precision * tpr) / (precision + tpr)  # F1分数

            print('\n'+'='*60)
            print(f"zscore门限: {z_thres:.1f}")
            print(f"FPR (假阳性率): {fpr:.4f}")
            print(f"TPR/Recall (召回率): {tpr:.4f}")
            print(f"Precision (精确率): {precision:.4f}")
            print(f"F1-Score: {f1_score:.4f}")
            cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis] * 100
            disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=["control", "dA-AL-II"])
            disp.plot(cmap=plt.cm.Blues, values_format='.2f', colorbar=False)
            plt.title('zscore threshold: %.1f' % z_thres)
            plt.savefig(os.path.join(fig_save_dir, 'zscore_threshold_%d_confusion_matrix.png' % z_thres), dpi=300)
            plt.close()

        # ROC曲线：zscore 作为连续分数，并标注 0.5-3.0（间隔0.5）的阈值点
        fpr, tpr, thresholds = roc_curve(y_true, feas)
        roc_auc = auc(fpr, tpr)
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, lw=2, label=f'ROC curve (AUC={roc_auc:.3f})')
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random Guess')
        for z_thres in np.arange(0.5, 2.6, 0.5):
            idx = np.argmin(np.abs(thresholds - z_thres))
            plt.scatter(fpr[idx], tpr[idx], s=50, zorder=5)
            plt.text(fpr[idx], tpr[idx], f' z={z_thres:.1f}', fontsize=9, va='bottom', ha='left')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Curve with z-score thresholds')
        plt.grid(alpha=0.3)
        plt.legend(loc='lower right')
        plt.savefig(os.path.join(fig_save_dir, 'roc_curve.png' % z_thres), dpi=300)

        # 按 0-1.0（间隔 0.05）比例混合正负样本，zscore 检测门限=z_thres_mix
        # 使用不同随机种子重复 100 次，绘制热力图
        n_base = int(min(len(control_feas)*0.8, len(dAAL_feas)*0.8))
        mix_props = np.arange(0, 1.01, 0.02)
        n_repeats = 100
        pred_props_matrix = np.zeros((n_repeats, len(mix_props)))

        for seed in range(n_repeats):
            np.random.seed(seed)
            for j, p in enumerate(mix_props):
                n_pos = int(round(p * n_base))
                n_neg = n_base - n_pos
                # 边界处理：比例为0或1时只取单类
                pos_sample = np.random.choice(dAAL_feas, size=n_pos, replace=False) if n_pos > 0 else np.array([])
                neg_sample = np.random.choice(control_feas, size=n_neg, replace=False) if n_neg > 0 else np.array([])
                mixed = np.concatenate([pos_sample, neg_sample])
                pred_props_matrix[seed, j] = np.mean(mixed >= z_thres_mix)

        # 每次重复的相关系数及整体相关系数
        corrs = np.array([np.corrcoef(mix_props, pred_props_matrix[i])[0, 1] for i in range(n_repeats)])
        overall_corr = np.corrcoef(np.tile(mix_props, n_repeats), pred_props_matrix.flatten())[0, 1]
        print('\n'+'='*60)
        print(f"zscore检测门限={z_thres_mix:.1f}，100 次随机抽样的相关系数: {corrs.mean():.4f} ± {corrs.std():.4f}")
        print(f"整体相关系数: {overall_corr:.4f}")

        # 绘制设置比例 vs 预测比例 的二维热力图
        pred_bins = np.arange(0, 1.01, 0.02)
        hist, xedges, yedges = np.histogram2d(
            np.tile(mix_props, n_repeats),   # x: 设置比例
            pred_props_matrix.flatten(),     # y: 预测比例
            bins=[mix_props, pred_bins]
        )

        # 将 count 转换为 ratio（占 100 次重复的比例）
        hist_ratio = hist / n_repeats

        # 横纵坐标 tick label 每隔 3 个显示一个
        x_labels = [f'{p:.2f}' if i % 3 == 0 else '' for i, p in enumerate(mix_props)]
        y_labels = [f'{p:.2f}' if i % 3 == 0 else '' for i, p in enumerate(pred_bins)]

        plt.figure(figsize=(10, 8))
        ax = sns.heatmap(hist_ratio.T, xticklabels=x_labels,
                        yticklabels=y_labels,
                        cmap='YlOrRd', cbar_kws={'label': 'Ratio'}, vmin=0, vmax=1)
        ax.invert_yaxis()
        # 叠加理想对角线 y=x（热力图单元中心对齐）
        # plt.plot([0.5, len(mix_props)-0.5], [0.5, len(pred_bins)-1.5], 'b--', lw=2, label='Ideal (y=x)')
        plt.xlabel('Set proportion of positive samples')
        plt.ylabel('Predicted positive proportion')
        plt.title(f'Set proportion vs Predicted proportion(z-score threshold={z_thres_mix:.1f}, overall corr={overall_corr:.3f})')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(fig_save_dir, 'hotmap.png' % z_thres), dpi=300)
        # plt.show()
        plt.close()
