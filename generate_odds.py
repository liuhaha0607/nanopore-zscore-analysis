#!/usr/bin/env python3
# ELIGOS2_complete.py - 完整版：Odds Ratio + ESB雷达图 + 质量值统计

import os
import pysam
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from collections import defaultdict


def get_mean_qscores(qvalues):
    """
    Convert qstring into a mean qscore
    """
    if len(qvalues) == 0: return 0.0
    err_probs = [10**(q / -10) for q in qvalues]
    mean_err = np.mean(err_probs)
    return -10 * np.log10(max(mean_err, 1e-4))

def calculate_error_rate_with_qstats(bam_file, ref_fasta, max_reads=10000, min_base_quality=20):
    """
    准确计算错误率，并添加质量值统计
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    ref_dict = pysam.FastaFile(ref_fasta)
    ref_name = ref_dict.references[0]
    ref_seq = ref_dict.fetch(ref_name)
    ref_len = len(ref_seq)
    
    # 初始化计数器
    total_reads = defaultdict(int)
    mismatch_reads = defaultdict(int)
    insertion_reads = defaultdict(int)
    deletion_reads = defaultdict(int)
    
    # 质量值统计
    pos_q = {}
    err_pos_q = {}
    
    base_counts = defaultdict(lambda: defaultdict(int))
    read_count = 0
    
    for read in bam.fetch(ref_name):
        if read_count >= max_reads:
            break
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        
        query_seq = read.query_sequence
        query_qual = read.query_qualities
        
        if query_seq is None or query_qual is None:
            continue
        
        aligned_pairs = read.get_aligned_pairs(matches_only=False)
        
        for query_pos, ref_pos in aligned_pairs:
            if ref_pos is None:
                continue
            
            if ref_pos < 0 or ref_pos >= ref_len:
                continue
            
            total_reads[ref_pos] += 1
            
            if query_pos is None:
                deletion_reads[ref_pos] += 1
                base_counts[ref_pos]['-'] += 1
            else:
                qual = query_qual[query_pos]
                
                if ref_pos not in pos_q:
                    pos_q[ref_pos] = [qual]
                else:
                    pos_q[ref_pos].append(qual)
                ref_base = ref_seq[ref_pos].upper()
                query_base = query_seq[query_pos].upper()
                base_counts[ref_pos][query_base] += 1
                
                is_mismatch = (ref_base != query_base)
                if is_mismatch:
                    mismatch_reads[ref_pos] += 1
                    if ref_pos not in err_pos_q:
                        err_pos_q[ref_pos] = [qual]
                    else:
                        err_pos_q[ref_pos].append(qual)
        
        read_count += 1
    
    bam.close()
    ref_dict.close()
    
    results = []
    for pos in range(ref_len):
        total = total_reads.get(pos, 0)
        mismatches = mismatch_reads.get(pos, 0)
        insertions = insertion_reads.get(pos, 0)
        deletions = deletion_reads.get(pos, 0)
        
        errors = mismatches + insertions + deletions
        error_rate = errors / total if total > 0 else 0
        
        if pos not in pos_q:
            q_mean = 0
        else:
            q_mean = get_mean_qscores(pos_q[pos])
        
        if pos not in err_pos_q:
            err_q_mean = 0
        else:
            err_q_mean = get_mean_qscores(err_pos_q[pos])

        results.append({
            'position': pos + 1,
            'ref_base': ref_seq[pos],
            'total_reads': total,
            'mismatches': mismatches,
            'insertions': insertions,
            'deletions': deletions,
            'total_errors': errors,
            'error_rate': error_rate,
            'mismatch_rate': mismatches / total if total > 0 else 0,
            'insertion_rate': insertions / total if total > 0 else 0,
            'deletion_rate': deletions / total if total > 0 else 0,
            'q_mean': q_mean,
            'error_q_mean': err_q_mean,
        })
    
    return pd.DataFrame(results)

def calculate_odds_ratio_with_qstats(test_bam, control_bam, ref_fasta, min_reads=50):
    """
    计算Odds Ratio，包含质量值统计
    """
    print("计算Test样本...")
    test_df = calculate_error_rate_with_qstats(test_bam, ref_fasta)
    
    print("计算Control样本...")
    control_df = calculate_error_rate_with_qstats(control_bam, ref_fasta)
    
    merge_cols = ['position', 'ref_base', 'total_reads', 'total_errors', 'error_rate',
                  'mismatches', 'deletions', 'insertions', 'q_mean', 'error_q_mean']
    
    merged = pd.merge(
        test_df[merge_cols],
        control_df[merge_cols],
        on=['position', 'ref_base'],
        suffixes=('_test', '_control')
    )
    
    merged = merged[
        (merged['total_reads_test'] >= min_reads) & 
        (merged['total_reads_control'] >= min_reads)
    ].copy()
    
    pseudocount = 0.5
    merged['odds_test'] = (merged['total_errors_test'] + pseudocount) / \
                          (merged['total_reads_test'] - merged['total_errors_test'] + pseudocount)
    merged['odds_control'] = (merged['total_errors_control'] + pseudocount) / \
                             (merged['total_reads_control'] - merged['total_errors_control'] + pseudocount)
    merged['odds_ratio'] = merged['odds_test'] / merged['odds_control']
    
    p_values = []
    for _, row in merged.iterrows():
        table = [
            [int(row['total_errors_test']), int(row['total_reads_test'] - row['total_errors_test'])],
            [int(row['total_errors_control']), int(row['total_reads_control'] - row['total_errors_control'])]
        ]
        try:
            _, pval = fisher_exact(table, alternative='greater')
        except:
            pval = 1.0
        p_values.append(pval)
    
    merged['p_value'] = p_values
    merged['-log10_p'] = -np.log10(np.array(p_values) + 1e-300)
    
    return merged

def plot_odds_ratio(df, output_file='odds_ratio.png', highlight_pos=None, title_suffix=""):
    """
    绘制Odds Ratio图（类似论文图4a,b）
    """
    fig, ax = plt.subplots(figsize=(14, 5))
    
    sig_df = df.query('p_value < 0.05 and error_rate_test > 0.1 and odds_ratio > 5')
    
    ax.scatter(df['position'], df['odds_ratio'], c='lightgray', s=10, alpha=0.5, label='All sites')
    
    if len(sig_df) > 0:
        ax.scatter(sig_df['position'], sig_df['odds_ratio'], 
                  c='red', s=15, alpha=0.8, label=f'Significant (p<0.05))')
    
    if highlight_pos:
        ax.axvline(x=highlight_pos, color='blue', linestyle='--', alpha=0.5)
    
    # ax.axhline(y=1, color='gray', linestyle='-', linewidth=0.5)
    # ax.axhline(y=2, color='orange', linestyle='--', alpha=0.3, label='OR=2')
    
    if len(sig_df) > 0:
        peak_pos = sig_df.loc[sig_df['odds_ratio'].idxmax(), 'position']
        zoom_start = max(0, int(peak_pos) - 30)
        zoom_end = min(int(df['position'].max()), int(peak_pos) + 30)
        
        axins = inset_axes(ax, width="35%", height="35%", loc='upper right')
        
        zoom_df = df[(df['position'] >= zoom_start) & (df['position'] <= zoom_end)]
        axins.scatter(zoom_df['position'], zoom_df['odds_ratio'], c='red', s=20)
        axins.plot(zoom_df['position'], zoom_df['odds_ratio'], 'r-', linewidth=1.5, alpha=0.7)
        axins.axhline(y=1, color='gray', linestyle='-', linewidth=0.5)
        axins.set_xlim(zoom_start, zoom_end)
        # axins.set_title(f'Zoom: {zoom_start}-{zoom_end} bp', fontsize=9)
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    
    ax.set_xlabel('Sequence Position (bp)', fontsize=12)
    ax.set_ylabel('Odds Ratio', fontsize=12)
    ax.set_title(f'Odds Ratio', fontsize=14)
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3, linestyle=':')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"已保存: {output_file}")
    plt.close()

def plot_esb_radar(df, center_pos, window=10, output_file='esb_radar.png', title_suffix=""):
    """
    绘制ESB雷达图（类似论文图4c,d）
    """
    window_df = df[(df['position'] >= center_pos - window) & 
                   (df['position'] <= center_pos + window)].copy()
    
    if len(window_df) < window * 2 + 1:
        print(f"窗口内数据不足: {len(window_df)} positions")
        return
    
    positions = window_df['position'].values
    angles = np.linspace(-1/4 * np.pi, 5/4 * np.pi,len(positions), endpoint=True)[::-1]
    angles = np.where(angles < 0, angles + 2 * np.pi, angles)
    
    test_errors = window_df['error_rate_test'].values.tolist()
    control_errors = window_df['error_rate_control'].values.tolist()
    
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='polar'))
    
    ax.plot(angles, test_errors, 'o-', linewidth=2, label='dA-AL-II', color='red')
    # ax.fill(angles, test_errors, alpha=0.25, color='red')
    
    ax.plot(angles, control_errors, 'o-', linewidth=2, label='control', color='black')
    # ax.fill(angles, control_errors, alpha=0.25, color='black')
    
    max_idx = np.argmax(test_errors)
    ax.plot(angles[max_idx], test_errors[max_idx], '*', markersize=15, color='red')
    
    labels = [f"{int(p)}:{row['ref_base']}" for p, (_, row) in zip(positions, window_df.iterrows())]
    ax.set_xticks(angles)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_rlabel_position(angles[0]*180/np.pi+5)
    ax.set_ylim(0, max(max(test_errors), max(control_errors)) * 1.2)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
    # ax.set_title(f'ESB around position {center_pos} bp{title_suffix}', y=1.08, fontsize=14)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"已保存: {output_file}")
    plt.close()

def plot_quality_comparison(df, output_file='quality_comparison.png'):
    """
    绘制质量值对比图
    """
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.plot(df['q_mean_test'], label='dA-AL-II', color='red')
    ax.plot(df['q_mean_control'], label='control', color='black')
    ax.set_xticks(np.arange(0, len(df), 40), np.arange(0, len(df), 40), rotation=90)
    ax.set_xlabel('position')
    ax.set_ylabel('Mean QValue')
    ax.grid()
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"已保存质量值对比图: {output_file}")
    # plt.show()
    plt.close()

# 主函数
if __name__ == "__main__":
    
    platform="ONT"
    root_dir = os.path.dirname(os.path.abspath(__file__))
    control_dir = os.path.join(root_dir, 'data', f'{platform}_control.bam')
    dAAL_dir = os.path.join(root_dir, 'data', f'{platform}_dAAL.bam')

    fig_save_dir = os.path.join(root_dir, 'fig', platform)
    data_save_dir = os.path.join(root_dir, 'data')
    os.makedirs(fig_save_dir, exist_ok=True)
    os.makedirs(data_save_dir, exist_ok=True)
    
    ref_fa = os.path.join(root_dir, 'data', 'ref.fasta')
    print("="*60)
    print("开始完整分析（Odds Ratio + ESB雷达图 + 质量值统计）")
    print("="*60)
    
    # 计算Odds Ratio
    results = calculate_odds_ratio_with_qstats(dAAL_dir, control_dir, ref_fa, min_reads=30)
    results.to_csv(os.path.join(data_save_dir, f"{platform}_odds_ratio.csv"), index=False)
    print(f"\n结果已保存到: {platform}_odds_ratio.csv")
    
    # 统计摘要
    print(f"\n{'='*60}")
    print("统计摘要")
    print(f"{'='*60}")
    print(f"总位点数: {len(results)}")
    print(f"\nTest样本:")
    print(f"  平均覆盖度: {results['total_reads_test'].mean():.1f}")
    print(f"  平均错误率: {results['error_rate_test'].mean():.4f}")
    print(f"  平均质量值: {results['q_mean_test'].mean():.2f}")
    print(f"  错误碱基平均质量: {results['error_q_mean_test'].mean():.2f}")
    
    print(f"\nControl样本:")
    print(f"  平均覆盖度: {results['total_reads_control'].mean():.1f}")
    print(f"  平均错误率: {results['error_rate_control'].mean():.4f}")
    print(f"  平均质量值: {results['q_mean_control'].mean():.2f}")
    print(f"  错误碱基平均质量: {results['error_q_mean_control'].mean():.2f}")
    
    # 显著位点
    sig_results = results[results['p_value'] < 0.05]
    print(f"\n显著位点 (p<0.05): {len(sig_results)} 个")
    
    if len(sig_results) > 0:
        print("\nTop 5 位点 (按Odds Ratio排序):")
        top5 = sig_results.nlargest(5, 'odds_ratio')[['position', 'ref_base', 'odds_ratio', 
                'p_value', 'error_rate_test', 'error_rate_control', 
                'q_mean_test', 'error_q_mean_test']]
        print(top5.to_string())
        
        peak_pos = int(sig_results.loc[sig_results['odds_ratio'].idxmax(), 'position'])
        
        # 生成所有图表
        print(f"\n{'='*60}")
        print("生成图表...")
        print(f"{'='*60}")
        
        # 1. Odds Ratio图
        plot_odds_ratio(results, os.path.join(fig_save_dir, f"{platform}_odds_ratio.png"), 
                       highlight_pos=peak_pos, title_suffix=" (with Indel Fix)")
        
        # 2. ESB雷达图
        plot_esb_radar(results, 506, window=5, 
                      output_file=os.path.join(fig_save_dir, f"{platform}_esb_radar.png"),
                      title_suffix=" (with Q-stats)")
        
        # 3. 质量值对比图
        plot_quality_comparison(results, os.path.join(fig_save_dir, f"{platform}_quality.png"))
        
        print(f"\n所有文件保存至: {fig_save_dir}")
        print(f"  - {platform}_odds_ratio.png")
        print(f"  - {platform}_esb_radar.png")  
        print(f"  - {platform}_quality.png")
        print(f"  - {platform}_odds_ratio.csv")
        
    else:
        print("未发现显著位点，仅生成基础图表")
        plot_odds_ratio(results, os.path.join(fig_save_dir, f"{platform}_odds_ratio.png"), 
                       title_suffix=" (No Significant Sites)")
    
    print(f"{'='*60}")
    print("分析完成!")
    print(f"{'='*60}")