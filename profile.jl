using Profile
using WaveToyDagger

# @time WaveToyDagger.main()

@profile WaveToyDagger.main()

open("profile.txt", "w") do fh
    return Profile.print(IOContext(fh, :displaysize => (24, 500)))
end
