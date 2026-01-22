# r-universe 配置文件

这个文件夹包含发布到 r-universe 所需的配置。

## 已完成配置

- 仓库: https://github.com/Zaoqu-Liu/zaoqu-liu.r-universe.dev
- r-universe App 已安装
- 包页面: https://zaoqu-liu.r-universe.dev

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
