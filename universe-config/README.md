# r-universe 配置文件

这个文件夹包含发布到 r-universe 所需的配置。

## 设置步骤

### 1. 创建 universe 仓库

在 GitHub 上创建新仓库: `https://github.com/Zaoqu-Liu/universe`

**重要**: 仓库名必须是 `universe`

### 2. 上传配置文件

将 `packages.json` 文件上传到该仓库的根目录。

### 3. 等待部署

r-universe 会在约 1 小时内自动检测并构建你的包。

### 4. 访问你的 r-universe

部署完成后访问: https://zaoqu-liu.r-universe.dev

## 用户安装方式

```r
install.packages("scMetaLink", repos = "https://zaoqu-liu.r-universe.dev")
```

## 添加更多包

如果以后有更多 R 包要发布，只需编辑 `packages.json`:

```json
[
  {
    "package": "scMetaLink",
    "url": "https://github.com/Zaoqu-Liu/scMetaLink"
  },
  {
    "package": "AnotherPackage",
    "url": "https://github.com/Zaoqu-Liu/AnotherPackage"
  }
]
```

## 自动更新

每次你 push 代码到 scMetaLink 仓库，r-universe 会自动重新构建并更新包。
