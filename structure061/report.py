import json
from pathlib import Path
p=Path(__file__).resolve().parent
rows=json.loads((p/'results/summary.json').read_text());old=json.loads((p.parent/'structure054/results/summary.json').read_text())
changes=[(a['id'],b['candidate'],a['candidate']) for a,b in zip(rows,old) if a['candidate']!=b['candidate']]
s='''# 055 名称归一化与复杂形貌入口

实验副本；未部署、未回填、未修改原始文件。LAM/LAMs归一为LAM，PtS/PtS1/PtS2归一为PtS；reference_key保存未归一的散射库候选，candidate_before_normalization保存最终候选原值，trials/ranking继续保留所有原库键。未改原始反射库。

通用形貌函数只接收计算所得拓扑和晶胞，不接受样本编号或文件名。三个阈值均有三个分别沿三轴贯穿的独立域时输出TriCylinder形貌候选；同时有两个孤立域且立方胞内位置相差(1/2,1/2,1/2)时输出TriCylinder+BCC。这些是形貌识别，不证明理想柱形、球形、空间群或平衡稳定性。

孤立域全部Euler数0才可输出Torus形貌；Euler数0与1共存则输出structure_name='Torus-containing mixed domains'，不赋予纯Torus最终相名。当前Torus实例属于后一种。Euler特征是离散拓扑证据，不是完整曲面同胚证明，仍需独立样本验证。

15个结构控制测试：三例目标的原始、平移、轴置换结果一致；FCC/A15负对照不误报为这三种形貌。两组别名映射测试通过。全部66例重新执行，包含历史已下架种子作为回归对照，不将它们重新加入图鉴。所有输入哈希一致，confirmed_label=null，database_written=false。

候选变化：
'''+str(changes)+'\n\n|实例|规范候选|原散射库键|结构描述|状态|\n|---|---|---|---|---|\n'
for r in rows:
 d=json.loads((p/'results'/(r['id']+'.json')).read_text());assert d['input_hashes_before']==d['input_hashes_after'];assert d['confirmed_label'] is None and not d['database_written']
 s+=f"|{r['id']}|{r['candidate'] or '—'}|{d.get('reference_key') or '—'}|{d.get('structure_name') or '—'}|{r['status']}|\n"
(p/'REPORT.md').write_text(s,encoding='utf-8');print(changes);print('Verified',len(rows))
