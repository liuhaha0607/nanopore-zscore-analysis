
import os
import numpy as np
import pickle
import scipy.signal as ss
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde

index_to_base = ['A', 'C', 'G', 'T']   
ref = 'AAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCCGCCTGGCATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGGTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAGCGGTTTGACTCACGGGGATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCAAAATCAACGGGACTTTCCAAAATGTCGTAACAACTCCGCCCCATTGACGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCTCTCTGGCTAACTAGAGAACCCACTGCTTACTGGCTTATCGAAATTAATACGACTCACTATAGGGAGACCCAAGCTGGCTAGCGTTTAAACTTAAGCTTGGTACCGAGCTCGGATCCACTAGTCCAGTGTGGTGGAATTCTGCAGATATCCAGCACAGTGGCGGCCGCGCGCGGGCCCAGCCGGCGCCCGTCCGCGCCGGGCCCGTTTAAACCCGCTGATCAGCCTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTCTGAGTAGGTGTCATTCTATTCTGGGGGGTGGGGTGGGGCAGGACAGCAAGGGGGAGGATTGGGAAGACAATAGCAGGCATGCTGGGGATGCGGTGGGCTCTATGGCTTCTGAGGCGGAAAGAACCAGCTGGGGCTCTAGGGGGTATCCCCACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACCCCAAAAAACTTGATTAGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTAATTCTGTGGAATGTGTGTCAGTTAGGGTGTGGAAAGTCCCCAGGCTCCCCAGCAGGCAGAAGTATGCAAAGCATGCATCTCAATTAGTCAGCAACCAGGTGTGGAAA'


if __name__ == '__main__':
    # 绘图参数
    colors = ['k', 'r']
    # colors = ['r', 'b']
    # colors = ['#20A23B', '#B81D2E']
    alpha_scaler = 10

    for platform in ['QT', 'ONT']:
        if platform=='ONT':
            start_position=11
        else:
            start_position=7
        root_dir = os.path.dirname(os.path.abspath(__file__))
        sig_dir_list = {
            'control': os.path.join(root_dir, 'data', f'{platform}_control_align_res.pkl'),
            'dA-AL-II': os.path.join(root_dir, 'data', f'{platform}_dAAL_align_res.pkl'),
        }

        fig_save_dir = os.path.join(root_dir, 'fig', platform)
        os.makedirs(fig_save_dir, exist_ok=True)

        # 加载数据
        all_sigs = {}
        resample = 6
        for chip_name, cache_save_path in sig_dir_list.items():
            with open(cache_save_path, 'rb') as fp:
                s = pickle.load(fp)
                cur_sigs = s['sigs']
            all_sigs[chip_name] = cur_sigs

        print('===绘图===')
        # 全局信号图
        fig, ax = plt.subplots(figsize=(20, 6))
        min_pos = 27
        max_pos = 1320
        for j, (chip_name, cur_sigs) in enumerate(all_sigs.items()):
            cur_sigs = np.array(cur_sigs)
            alpha = 1 / len(cur_sigs) * alpha_scaler  # 线条透明度， 根据read数量自适应设定
            for i in range(len(cur_sigs)):
                cur_sig = ss.medfilt(cur_sigs[i][(min_pos-start_position)*resample: (max_pos-start_position)*resample], kernel_size=5)
                if i == 0:
                    ax.plot(np.arange(min_pos*resample, min_pos*resample+1)/resample, [-10], color=colors[j], linewidth=1, alpha=1.0, label=chip_name)
                    ax.plot(np.arange(min_pos*resample, max_pos*resample)/resample, cur_sig, color=colors[j], linewidth=1, alpha=2*alpha)
                else:
                    ax.plot(np.arange(min_pos*resample, max_pos*resample)/resample, cur_sig, color=colors[j], linewidth=1, alpha=alpha)
        ax.set_ylim(-3.0, 3.0)
        ax.set_ylabel('signal')
        ax.set_xlabel('position')
    
        ax.set_xticks(np.arange(min_pos, max_pos, 40)-1/resample, np.arange(min_pos, max_pos, 40)+1, rotation=0) # 此处+1是为了从1开始计数
        plt.vlines(np.array([111, 121, 501, 511]), -3, 3, color='k', linestyle='--')
        plt.savefig(os.path.join(fig_save_dir, 'sig_compare_all.png'), dpi=300)

        # 非加合段信号图
        fig, ax = plt.subplots(figsize=(10, 4))
        min_pos = 111
        max_pos = 121
        for j, (chip_name, cur_sigs) in enumerate(all_sigs.items()):
            cur_sigs = np.array(cur_sigs)
            alpha = 1 / len(cur_sigs) * alpha_scaler
            for i in range(len(cur_sigs)):
                cur_sig = ss.medfilt(cur_sigs[i][(min_pos-start_position)*resample: (max_pos-start_position)*resample], kernel_size=5)
                if i == 0:
                    ax.plot(np.arange(min_pos*resample, min_pos*resample+1)/resample, [-10], color=colors[j], linewidth=1, alpha=1.0, label=chip_name)
                    ax.plot(np.arange(min_pos*resample, max_pos*resample)/resample, cur_sig, color=colors[j], linewidth=1, alpha=2*alpha)
                else:
                    ax.plot(np.arange(min_pos*resample, max_pos*resample)/resample, cur_sig, color=colors[j], linewidth=1, alpha=alpha)
        
        text_base = []
        for i in np.arange(min_pos, max_pos):
            text_base.append(ref[i])

        ax.legend()
        ax.set_ylim(-3.0, 3.0)
        ax.set_ylabel('signal')
        ax.set_xlabel('position')
    
        ax.set_xticks(np.arange(min_pos, max_pos)-1/resample, np.arange(min_pos, max_pos)+1, rotation=0)
        # 主网格线（对应主刻度，每5个单位）
        ax.grid(which='major', linestyle='-', linewidth=0.8, color='gray', alpha=0.7)
        # 次网格线（对应次刻度，每1个单位）
        for loc, base in zip(np.arange(min_pos, max_pos), text_base):
            ax.text(loc+0.2, 2.5, base, color='r' if int(loc)==505 else 'k')
        plt.savefig(os.path.join(fig_save_dir, 'sig_compare_font.png'), dpi=300)

        # 非加合段分布
        fig, axes = plt.subplots(1, max_pos-min_pos, figsize=(14, 6), sharey=True)
        for j, ax in enumerate(axes):
            # 提取当前位置 j 的所有信号
            pos_idx = j + min_pos
            
            for i, (chip_name, cur_sigs) in enumerate(all_sigs.items()):
                cur_sigs = np.array(cur_sigs)  # shape: (n_replicates, n_samples, signal_length)
                # 提取当前位置的信号片段
                start_idx = (pos_idx - start_position) * resample
                end_idx = (pos_idx + 1 - start_position) * resample
                cur_sig = cur_sigs[:, start_idx:end_idx].flatten()
                if len(cur_sig) > 0:
                    # KDE直接对原始数据
                    kde = gaussian_kde(cur_sig, bw_method=0.3)
                    # 创建评估点
                    y_range = np.linspace(cur_sig.min() - 0.5, cur_sig.max() + 0.5, 500)
                    density = kde(y_range)
                    ax.plot(density, y_range, color=colors[i], linewidth=1.5, label=chip_name)
            # 设置样式
            ax.axvline(x=0, color='gray', linewidth=0.8)
            ax.set_ylim(-4, 4)

            if pos_idx+1 == 506:
                ax.set_title(f'{pos_idx+1}:{ref[pos_idx]}', color='r')
            else:
                ax.set_title(f'{pos_idx+1}:{ref[pos_idx]}', color='k')
            
            # 只在第一个子图显示y轴标签和图例
            if j == 0:
                ax.set_ylabel('Signal', fontsize=12)
                ax.legend(loc='upper right')
            else:
                ax.set_yticklabels([])
            ax.set_xticklabels([])
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.05)
        plt.savefig(os.path.join(fig_save_dir, 'sig_distribution_font.png'), dpi=300)

        # 加和段信号图
        fig, ax = plt.subplots(figsize=(10, 4))
        min_pos = 501
        max_pos = 511
        for j, (chip_name, cur_sigs) in enumerate(all_sigs.items()):
            cur_sigs = np.array(cur_sigs)
            alpha = 1 / len(cur_sigs) * alpha_scaler
            for i in range(len(cur_sigs)):
                cur_sig = ss.medfilt(cur_sigs[i][(min_pos-start_position)*resample: (max_pos-start_position)*resample], kernel_size=5)
                if i == 0:
                    ax.plot(np.arange(min_pos*resample, min_pos*resample+1)/resample, [-10], color=colors[j], linewidth=1, alpha=1.0, label=chip_name)
                    ax.plot(np.arange(min_pos*resample, max_pos*resample)/resample, cur_sig, color=colors[j], linewidth=1, alpha=2*alpha)
                else:
                    ax.plot(np.arange(min_pos*resample, max_pos*resample)/resample, cur_sig, color=colors[j], linewidth=1, alpha=alpha)
        
        text_base = []
        for i in np.arange(min_pos, max_pos):
            text_base.append(ref[i])

        ax.legend()
        ax.set_ylim(-3.0, 3.0)
        ax.set_ylabel('signal')
        ax.set_xlabel('position')
    
        ax.set_xticks(np.arange(min_pos, max_pos)-1/resample, np.arange(min_pos, max_pos)+1, rotation=0)
        # 主网格线（对应主刻度，每5个单位）
        ax.grid(which='major', linestyle='-', linewidth=0.8, color='gray', alpha=0.7)
        # 次网格线（对应次刻度，每1个单位）
        for loc, base in zip(np.arange(min_pos, max_pos), text_base):
            ax.text(loc+0.2, 2.5, base, color='r' if int(loc)==505 else 'k')
        plt.savefig(os.path.join(fig_save_dir, 'sig_compare_dAAL.png'), dpi=300)

        # 加合段分布
        fig, axes = plt.subplots(1, max_pos-min_pos, figsize=(14, 6), sharey=True)
        for j, ax in enumerate(axes):
            # 提取当前位置 j 的所有信号
            pos_idx = j + min_pos
            
            for i, (chip_name, cur_sigs) in enumerate(all_sigs.items()):
                cur_sigs = np.array(cur_sigs)  # shape: (n_replicates, n_samples, signal_length)
                # 提取当前位置的信号片段
                start_idx = (pos_idx - start_position) * resample
                end_idx = (pos_idx + 1 - start_position) * resample
                cur_sig = cur_sigs[:, start_idx:end_idx].flatten()
                if len(cur_sig) > 0:
                    # KDE直接对原始数据
                    kde = gaussian_kde(cur_sig, bw_method=0.3)
                    # 创建评估点
                    y_range = np.linspace(cur_sig.min() - 0.5, cur_sig.max() + 0.5, 500)
                    density = kde(y_range)
                    ax.plot(density, y_range, color=colors[i], linewidth=1.5, label=chip_name)
            # 设置样式
            ax.axvline(x=0, color='gray', linewidth=0.8)
            ax.set_ylim(-4, 4)

            if pos_idx+1 == 506:
                ax.set_title(f'{pos_idx+1}:{ref[pos_idx]}', color='r')
            else:
                ax.set_title(f'{pos_idx+1}:{ref[pos_idx]}', color='k')
            
            # 只在第一个子图显示y轴标签和图例
            if j == 0:
                ax.set_ylabel('Signal', fontsize=12)
                ax.legend(loc='upper right')
            else:
                ax.set_yticklabels([])
            ax.set_xticklabels([])
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.05)
        plt.savefig(os.path.join(fig_save_dir, 'sig_distribution_dAAL.png'), dpi=300)
            
        plt.show()
        plt.close()