# quals = "@@CFFFDDHHGHHJIJJIJIJJJIIJJJJGIJJJFGIGHHHFFFFEEEE@BDDDDDDDDDBDBDEDDDCACCDACCDDDDBDDDDDDDDDDBDC??CDDDDDBBBB@<CDDDDDDCCDACCC>ACCDDDDDCDAAACB8<55@BDD7<?B?"
# quals = "CCCFFFFFHHHHHJJJEFHGGIJIJGHJJJJJIJJJJJJIIJJJJJJJJJHHHHFFFFEECCCCDDCDDDCDDDBDDDDDDDCCDDCDBDDDDDBBB0<@CDEEDBBB@@0<5?>@C1?AA8<ABCDB?A?CC>:BCDD+4<13>A:3:>>"
# quals = "@@CDFFFFHGGGHIIGIIICADHBHIIIIIDHGHGIIIIIGIIIIIIIIIGIIIIIIII8BEBDDDFFFDDDB;BDDDCCDCDCDDDCCDCDDD@BB>?CA>C@BCDDBD8?>CCD<<AAB>>9>0<@5902:144>>@9<&)8&&9(>(("
# quals = "CCCFFFFFHHHGHJJJJHJIJJIJJIJJJIIJJJIGIJJJJJJJJJJJJJJJJJJJIJIJF=;=;;.;?=',(,(53(5(,((,(,5<C(>34(3&))55??B&))&)0(+((((0))00((((+((+&&&)&+(+((+((4+4((((()&"
# quals = "@@@FFFDFHHHFHGBGIIJJDIIIHJJJJJ@FHIGIJIJGIJJIIHGCHFFFFDD5&)&007((+((((((((+(+><)(29(+4(((+4?(&)))0(((((44(()5&&&&&&((((++4(+(+44:9(&++3((+&)&&&++&&&)(()"
# quals = "@@CFFFFDHDFFDHGIIIHGJJGHIIJFGGGGEHHIIGIFFEHFFEFDDDD8?:A@CDCDC8ADBCDDDCDDCDDDDBDCDDDCCCD?BB<:BC@C<>9809BBBDBCD@B9<@BDBD((>:99B<&0<9>9&88<B3)09>7)05)(42<"
# quals = "@@@FFFFDFHHHHGIIIIGGGHIGIJIIGEHHBEGI=GHEHHIEAEEGIJJHJ@EEAEFDA?B7>ACDACBBDDB@BB>BDBDDDDBDDCB99ACB?ACCCA>@A@BBBACCDCC058A(:::@C@DC@CC@:4>@C>>AAABB:4:@8?9"
# quals = "CCCFFFFFHHHHHIIJJJJJJIIJJJJJJJJJJJJGIGJJJHFFDEEEEEDDDDEDDDDCDDDBEEDDDDDDDDDDDDDDDCDDDDDDBBDDDDDBEDCD<BDDDDDDDCCDD@BDDDDCCBD>BB99A+:@DDCCD32>><8?B>?AD34"
# quals = "CCCFFFFFHHHHHJJJJJJJDHIJIJJJIJJGJJIJIJHFDDDDDDDEDDDDDDDDDDDEDDDDDDDDDDCDDDDDDDDDBDDDDDDDDBDDBD>CDCDDDDDDDDCC?CDDCDDCCABD@?BCDDDDEDDDDDDDD3:@AA?BB>05<<A"
# quals = "CCCFFFFFHHGHHJIJJJJJJJJJJIJJJJJIIIIIJJJJJIJJJJIIIIIHIIHHHHHFFFFFFCC>CEEDDCDDDCDDBDDDDDBDBDDDDDDDDDDDDDDDDCDE>>ACBACCB?CCACDDBDD?B5>9&0@:ACACC4<:A<+4:@A"
# quals = "CCCFFFFFHGGHHJJJHIJJJJGIIJJIIJIJHFFFFDDDDDDDDDDDCDDDDDDDEDDDDDDBDDDDDDCDDDDDDBBDDDBDBDDDDD5>BCCEDD@CDBD<7.?CC@CBDDBB>@BBDBCACC:+:8:>?BB3>?72<<4+4>825&&"
# quals = "CCCFFFFFHHGHHJIJJIJJIJJJJJJEHGGIGHFFDDDDDDDDDEDDDDDDDDDDDDDDDBBBDDDDDDDDDDD@B@BDDDBBB>B???@>?B@9ACCDDDDD?CCACCD8?B?ABDDBDDBD>@@DCC>:>9CCCC?ACC2?CA8A@(4"
quals = "@@<DDDD?FHFBBHIIIIEGI<FHIGGIDHGIIIIIFEII>FHIIFED@CBB=?ACCC?CBBBB;8<B<A8C@+:(8A@CCDCCC::::??9A@C:@C@CCCCCCCCBB@>BBBBCCB?<59>@9BCC(:@4>@A>CCC8<.9<<B<:(:>"

# Turn quality string into ascii values.
ord <- function(x) {
  s <- x[1]
  
  intbits <- .Machine$sizeof.long * 4
  bits <- raw(nchar(s) * intbits)
  
  r <- charToRaw(x)
  idx <- if (.Platform$endian == "little") 1:8 else 25:32
  for (i in 1:nchar(s)) {
    bits[idx + ((i - 1) * intbits)] <- rawToBits(r[i])
  }
  packBits(bits, 'integer')
}

# Normalize the quality values.
quals = ord(quals) - 33
x = seq(1,length(quals))

# Store the 6 types of regression used.
colors = rainbow(6)
degrees = c(1,3,5,7,9)

pdf("example_regression.pdf", width=12, height=7.5)
plot(x,quals, xlab = "read position", ylab = "quality value")
curr_color <- 1

# Add regression (0).
abline(h = mean(quals), col = colors[[1]], lwd = 2)
for (degree in degrees) {
  fit = lm(quals ~ poly(x, degree))
  lines(fit$fitted.values, col = colors[[curr_color + 1]], lty = curr_color + 1, lwd = 2)
  curr_color <- curr_color + 1
}
legend("bottomleft",
       c("regression (0)","regression (1)", "regression (3)", "regression (5)", "regression (7)", "regression (9)"),
       lty=seq(1,6),
       lwd=c(3.5,3.5), col=colors)

dev.off()
