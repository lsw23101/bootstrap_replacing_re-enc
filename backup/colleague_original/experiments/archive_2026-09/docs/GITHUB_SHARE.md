# Sharing this repository on GitHub

From the repository root:

```bash
git init
git add .
git commit -m "Initial CKKS unstable-F experiment archive"
```

If GitHub CLI is installed and authenticated, a private repository can be created with:

```bash
gh repo create ckks-unstable-f-control --private --source=. --remote=origin --push
```

To share it with a colleague, add them as a collaborator in the GitHub repository settings. Use a private repository until any unpublished research details are ready for public release.

