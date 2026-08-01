# GitHub 网页上传操作指南

## 第一步：下载并解压

1. 下载 `nanopore-zscore-analysis-ready-to-upload.zip`。
2. 在电脑上解压。
3. 打开解压后的 `nanopore-zscore-analysis` 文件夹。
4. 确认里面能看到 4 个 `.py` 文件、`README.md`、`requirements.txt` 等文件。

## 第二步：注册或登录 GitHub

1. 打开 https://github.com/
2. 没有账号时点击 **Sign up** 注册。
3. 有账号时点击 **Sign in** 登录。
4. 建议在账户设置中验证邮箱。

## 第三步：创建公开仓库

1. 登录后，点击右上角的 **+**。
2. 点击 **New repository**。
3. Repository name 填：`nanopore-zscore-analysis`
4. Description 可填：
   `Analysis code for nanopore signal comparison and Z-score-based detection.`
5. 选择 **Public**。
6. 不要勾选 **Add a README file**，因为压缩包内已经有 README。
7. 暂时不要选择 `.gitignore` 或 License，压缩包内已有 `.gitignore`。
8. 点击 **Create repository**。

## 第四步：上传文件

1. 在新仓库页面点击 **uploading an existing file**；若没有看到，点击
   **Add file → Upload files**。
2. 打开电脑上解压后的项目文件夹。
3. 选择文件夹内的全部内容，而不是选择外层文件夹本身。
4. 将全部文件和 `data`、`fig` 文件夹拖入 GitHub 上传区域。
5. 等待文件列表全部显示。
6. 在 Commit changes 区域：
   - 标题填写：`Initial release of analysis code`
   - 说明可留空。
7. 选择 **Commit directly to the main branch**。
8. 点击 **Commit changes**。

## 第五步：检查仓库

仓库首页应看到：

- 4 个 Python 脚本；
- README.md；
- requirements.txt；
- environment.yml；
- data 文件夹；
- fig 文件夹。

向下滚动时，GitHub 会自动显示 README 内容。

## 第六步：补充文章信息

打开 `README.md`，点击铅笔图标编辑，把以下占位符替换掉：

- `[ADD ACCESSION NUMBER]`
- `[ADD RUN ACCESSION NUMBERS]`
- `[ADD ZENODO DOI]`（获得 DOI 后再填）

编辑完后点击 **Commit changes**。

## 第七步：添加许可证

1. 在仓库首页点击 **Add file → Create new file**。
2. 文件名输入 `LICENSE`。
3. 点击右侧 **Choose a license template**。
4. 与共同作者和单位确认后，可选择 **MIT License**。
5. 填写年份和版权所有者。
6. 点击 **Review and submit → Commit changes**。

## 第八步：发布固定版本

1. 在仓库首页右侧点击 **Releases**。
2. 点击 **Draft a new release**。
3. 点击 **Choose a tag**，输入 `v1.0.0`。
4. 点击 **Create new tag: v1.0.0 on publish**。
5. Release title 填：`Analysis code used in the manuscript`
6. Release notes 可填：
   `Version 1.0.0 of the analysis scripts used to generate the results and figures reported in the manuscript.`
7. 点击 **Publish release**。

然后可连接 Zenodo，为该固定版本生成 DOI。
