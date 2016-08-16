# Make sure your are up-to-date
git pull https://github.com/christopherwharrop/pasilla.git

# Look to see what files you've modified and have staged for committing
git status

# Add files to be committed
git add <list of files you want to commit>

# Do the commit
git commit -m "Description of the change"  # For a short message
git commit                                 # When a longer message is required

# Push you commit to the repository
git push origin master
