#!/usr/bin/env bash
#
# Push this repository to GitHub over HTTPS.
#
# GitHub has not accepted account passwords for git operations since August
# 2021.  At the "Password for ..." prompt, paste a Personal Access Token:
#
#     github.com -> Settings -> Developer settings -> Personal access tokens
#     -> Tokens (classic) -> Generate new token -> tick the "repo" scope
#
# Create the repository on github.com first -- empty, with no README,
# .gitignore or licence, or the first push will be rejected as diverged --
# then run:
#
#     ./push_to_github.sh YOUR-GITHUB-USERNAME [repo-name]
#
set -euo pipefail

user="${1:-}"
repo="${2:-dragonrna-design}"

if [ -z "$user" ]; then
    echo "usage: $0 YOUR-GITHUB-USERNAME [repo-name]" >&2
    echo "       (repo-name defaults to dragonrna-design)" >&2
    exit 1
fi

# Putting the username in the URL means git only asks for the token.
url="https://${user}@github.com/${user}/${repo}.git"

if git remote get-url origin >/dev/null 2>&1; then
    git remote set-url origin "$url"
else
    git remote add origin "$url"
fi

branch="$(git symbolic-ref --short HEAD)"

echo "Pushing branch '${branch}' to ${url}"
echo
echo "When prompted for a password, paste your Personal Access Token."
echo "It will not echo to the screen. Your account password will not work."
echo

git push -u origin "$branch"

echo
echo "Done: https://github.com/${user}/${repo}"
