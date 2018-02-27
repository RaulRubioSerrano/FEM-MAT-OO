%==================================================================
%                        General Data File
% Title: strip_hole
% Units: SI
% Dimensions: 2D
% Type of problem: Plane_Stress
% Type of Phisics: HYPERELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'TRIANGLE';
'SI';
'2D';
'Plane_Stress';
'HYPERELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z

coord = [
1          6.5          6.5            0
2       6.0357          6.5            0
3          6.5       6.0357            0
4       5.8314       6.0329            0
5          6.5       5.5714            0
6       5.5714          6.5            0
7       6.0838       5.5816            0
8        5.283       5.9888            0
9       5.5143       5.5912            0
10       5.1071          6.5            0
11          6.5       5.1071            0
12       5.7455       5.1936            0
13       5.0543       5.5897            0
14        4.823       5.9873            0
15       5.2855       5.1921            0
16       6.0364       4.7555            0
17       4.6429          6.5            0
18          6.5       4.6429            0
19       5.5168       4.7944            0
20       4.5943       5.5882            0
21       6.2081       4.3983            0
22       4.8256       5.1906            0
23        4.363       5.9859            0
24       5.7481       4.3968            0
25       5.0568        4.793            0
26       4.1786          6.5            0
27          6.5       4.1786            0
28       5.2881       4.3953            0
29       6.0377       4.0476            0
30       4.3656       5.1891            0
31       4.1343       5.5867            0
32       4.5968       4.7915            0
33        3.903       5.9844            0
34       5.5194       3.9977            0
35       4.8281       4.3939            0
36       3.7143          6.5            0
37          6.5       3.7143            0
38       5.0594       3.9962            0
39       3.9056       5.1876            0
40       4.1368         4.79            0
41       3.6743       5.5853            0
42       5.7507       3.6001            0
43       4.3681       4.3924            0
44        3.443       5.9829            0
45       5.2907       3.5986            0
46       4.5994       3.9947            0
47         3.25          6.5            0
48          6.5         3.25            0
49       6.0876       3.2257            0
50       3.6768       4.7885            0
51       3.4456       5.1862            0
52       3.9081       4.3909            0
53       4.8307       3.5971            0
54       3.2143       5.5838            0
55        5.522        3.201            0
56        4.098       3.9756            0
57        2.983       5.9814            0
58        5.062       3.1995            0
59       4.3707       3.5956            0
60       3.2168        4.787            0
61       2.7857          6.5            0
62          6.5       2.7857            0
63       3.4532       4.3454            0
64       2.9856       5.1847            0
65        5.836       2.7969            0
66       3.6794       3.9918            0
67        4.602        3.198            0
68       2.7543       5.5823            0
69       5.2932       2.8019            0
70        2.523       5.9799            0
71       4.1258       3.2257            0
72       4.8332       2.8004            0
73       2.9881       4.3879            0
74       2.7568       4.7856            0
75       3.6623       3.4881            0
76       6.0889       2.3631            0
77       2.3214          6.5            0
78          6.5       2.3214            0
79       2.5256       5.1832            0
80       5.5245       2.4042            0
81       4.3314       2.8346            0
82         3.25       3.7261            0
83       2.2943       5.5808            0
84       5.0645       2.4027            0
85        3.375       3.4665            0
86          3.5         3.25            0
87        2.063       5.9785            0
88       2.7909       3.9918            0
89       3.6624       3.0119            0
90       2.5281       4.3865            0
91       4.6045       2.4013            0
92        3.125       3.4665            0
93       2.2968       4.7841            0
94       5.7558       2.0066            0
95       2.0656       5.1817            0
96       1.8571          6.5            0
97          6.5       1.8571            0
98       5.2958       2.0051            0
99        3.375       3.0335            0
100       4.1445       2.3998            0
101       2.8377       3.4881            0
102       1.8343       5.5793            0
103            3         3.25            0
104       4.8358       2.0036            0
105        3.125       3.0335            0
106       2.3506       3.9704            0
107       6.0082       1.6148            0
108       2.0681        4.385            0
109        1.603        5.977            0
110         3.25       2.7739            0
111       1.8368       4.7826            0
112       4.3758       2.0022            0
113       5.5271       1.6075            0
114       3.6088       2.4304            0
115       2.8376       3.0119            0
116       1.6056       5.1802            0
117       5.0671        1.606            0
118       1.3929          6.5            0
119          6.5       1.3929            0
120       3.9158       2.0007            0
121       1.3743       5.5779            0
122       4.6071       1.6045            0
123       1.8394       3.9859            0
124       2.0707       3.5883            0
125        3.215       2.3193            0
126       2.3658       3.1437            0
127       1.6081       4.3835            0
128       5.7584       1.2099            0
129        1.143       5.9755            0
130       1.3769       4.7811            0
131       5.2984       1.2084            0
132       4.1471        1.603            0
133       3.4558       1.9992            0
134       2.4692       2.7473            0
135       1.1456       5.1788            0
136       4.8384       1.2069            0
137       2.7645       2.3953            0
138          6.5      0.92857            0
139      0.92857          6.5            0
140       6.0915      0.86991            0
141       3.6871       1.6016            0
142      0.91429       5.5764            0
143       1.6107       3.5868            0
144       4.3784       1.2054            0
145       1.3794       3.9844            0
146       2.9958       1.9977            0
147        1.842       3.1892            0
148       1.1481        4.382            0
149       5.5296      0.81075            0
150       2.0733       2.7915            0
151      0.68301        5.974            0
152      0.91685       4.7796            0
153       5.0696      0.80927            0
154       2.3045       2.3939            0
155       3.9184       1.2039            0
156       3.2271       1.6001            0
157       4.6096      0.80779            0
158       2.5358       1.9962            0
159      0.64798       5.1658            0
160          6.5      0.46429            0
161      0.46429          6.5            0
162       1.1507       3.5853            0
163        1.382       3.1877            0
164       3.4584       1.2025            0
165       5.8382       0.4261            0
166      0.45429       5.5749            0
167      0.91942       3.9829            0
168       1.6133         2.79            0
169       4.1496      0.80631            0
170       2.7671       1.5986            0
171       5.3009      0.41164            0
172       1.8446       2.3924            0
173      0.45686       4.7782            0
174      0.55724       4.3894            0
175       4.8409      0.41016            0
176       2.0758       1.9948            0
177       3.6896      0.80483            0
178       2.9984        1.201            0
179       4.3809      0.40868            0
180       2.3071       1.5971            0
181      0.92198       3.1862            0
182            0          6.5            0
183          6.5            0            0
184       1.1533       2.7885            0
185       6.0357            0            0
186            0       6.0357            0
187      0.65225       3.6058            0
188       1.3846       2.3909            0
189       5.5714            0            0
190            0       5.5714            0
191       3.2296      0.80335            0
192      0.42578       3.9742            0
193       3.9209       0.4072            0
194       2.5384       1.1995            0
195       1.6159       1.9933            0
196            0       5.1071            0
197       5.1071            0            0
198            0       4.6429            0
199       4.6429            0            0
200       1.8471       1.5957            0
201       3.4609      0.40572            0
202       2.7696      0.80188            0
203      0.46198       3.1847            0
204      0.69326       2.7871            0
205            0       4.1786            0
206       4.1786            0            0
207       2.0784        1.198            0
208      0.92454       2.3894            0
209       1.1559       1.9918            0
210       3.0009      0.40424            0
211            0       3.7143            0
212       3.7143            0            0
213       2.3097      0.80043            0
214       1.3872       1.5942            0
215       1.6184       1.1965            0
216         3.25            0            0
217            0         3.25            0
218        2.541      0.40282            0
219      0.46454        2.388            0
220       1.8497      0.79891            0
221      0.65659       1.9684            0
222      0.92726       1.5927            0
223            0       2.7857            0
224       2.7857            0            0
225       1.1585       1.1951            0
226        2.081       0.4013            0
227       1.3896      0.79744            0
228            0       2.3214            0
229       2.3214            0            0
230      0.42969       1.5991            0
231        1.621      0.39976            0
232      0.56462       1.1845            0
233      0.92963      0.79614            0
234       1.8571            0            0
235            0       1.8571            0
236       1.1607      0.39837            0
237            0       1.3929            0
238       1.3929            0            0
239      0.46963      0.79489            0
240      0.65881      0.40895            0
241      0.92857            0            0
242            0      0.92857            0
243            0      0.46429            0
244      0.46429            0            0
245            0            0            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Material

connec = [
1 86 99 89 0
2 105 103 115 0
3 92 85 82 0
4 99 105 110 0
5 103 92 101 0
6 85 86 75 0
7 89 99 110 0
8 82 85 75 0
9 101 92 82 0
10 115 103 101 0
11 110 105 115 0
12 75 86 89 0
13 182 186 161 0
14 241 238 236 0
15 236 238 231 0
16 231 238 234 0
17 236 231 227 0
18 227 231 220 0
19 236 227 233 0
20 220 231 226 0
21 227 220 215 0
22 236 233 240 0
23 220 226 213 0
24 215 220 207 0
25 240 233 239 0
26 213 226 218 0
27 215 207 200 0
28 239 233 232 0
29 218 226 229 0
30 218 229 224 0
31 200 207 180 0
32 232 233 225 0
33 239 232 242 0
34 180 207 194 0
35 200 180 176 0
36 225 233 227 0
37 225 227 215 0
38 225 215 214 0
39 214 215 200 0
40 225 214 222 0
41 214 200 195 0
42 195 200 176 0
43 214 195 209 0
44 195 176 172 0
45 194 207 213 0
46 213 207 220 0
47 194 213 202 0
48 180 194 170 0
49 176 180 158 0
50 209 195 188 0
51 170 194 178 0
52 176 158 154 0
53 188 195 172 0
54 178 194 202 0
55 154 158 137 0
56 188 172 168 0
57 178 202 191 0
58 137 158 146 0
59 168 172 150 0
60 191 202 210 0
61 146 158 170 0
62 150 172 154 0
63 154 172 176 0
64 150 154 134 0
65 210 202 218 0
66 218 202 213 0
67 210 218 224 0
68 210 224 216 0
69 170 158 180 0
70 146 170 156 0
71 156 170 178 0
72 146 156 133 0
73 156 178 164 0
74 133 156 141 0
75 164 178 191 0
76 156 164 141 0
77 141 164 155 0
78 155 164 177 0
79 141 155 132 0
80 177 164 191 0
81 155 177 169 0
82 132 155 144 0
83 177 191 201 0
84 169 177 193 0
85 144 155 169 0
86 201 191 210 0
87 193 177 201 0
88 169 193 179 0
89 201 210 216 0
90 201 216 212 0
91 179 193 206 0
92 179 206 199 0
93 169 179 157 0
94 157 179 175 0
95 169 157 144 0
96 175 179 199 0
97 175 199 197 0
98 157 175 153 0
99 144 157 136 0
100 153 175 171 0
101 157 153 136 0
102 171 175 197 0
103 171 197 189 0
104 136 153 131 0
105 131 153 149 0
106 136 131 117 0
107 149 153 171 0
108 131 149 128 0
109 117 131 113 0
110 136 117 122 0
111 128 149 140 0
112 131 128 113 0
113 122 117 104 0
114 140 149 165 0
115 113 128 107 0
116 104 117 98 0
117 122 104 112 0
118 165 149 171 0
119 107 128 119 0
120 98 117 113 0
121 112 104 91 0
122 165 171 189 0
123 165 189 185 0
124 98 113 94 0
125 91 104 84 0
126 94 113 107 0
127 98 94 80 0
128 84 104 98 0
129 91 84 72 0
130 80 94 76 0
131 98 80 84 0
132 72 84 69 0
133 76 94 97 0
134 84 80 69 0
135 72 69 58 0
136 69 80 65 0
137 58 69 55 0
138 72 58 67 0
139 69 65 55 0
140 58 55 45 0
141 67 58 53 0
142 55 65 49 0
143 58 45 53 0
144 49 65 62 0
145 53 45 38 0
146 38 45 34 0
147 53 38 46 0
148 34 45 42 0
149 38 34 28 0
150 46 38 35 0
151 42 45 55 0
152 28 34 24 0
153 35 38 28 0
154 24 34 29 0
155 28 24 19 0
156 35 28 25 0
157 29 34 42 0
158 19 24 16 0
159 25 28 19 0
160 29 42 37 0
161 16 24 21 0
162 25 19 15 0
163 21 24 29 0
164 16 21 18 0
165 18 21 27 0
166 27 21 29 0
167 15 19 12 0
168 12 19 16 0
169 15 12 9 0
170 12 16 11 0
171 9 12 7 0
172 7 12 11 0
173 9 7 4 0
174 4 7 3 0
175 9 4 8 0
176 8 4 6 0
177 9 8 13 0
178 13 8 14 0
179 9 13 15 0
180 14 8 10 0
181 13 14 20 0
182 15 13 22 0
183 20 14 23 0
184 22 13 20 0
185 15 22 25 0
186 23 14 17 0
187 25 22 32 0
188 32 22 30 0
189 25 32 35 0
190 30 22 20 0
191 32 30 40 0
192 35 32 43 0
193 30 20 31 0
194 40 30 39 0
195 43 32 40 0
196 31 20 23 0
197 39 30 31 0
198 43 40 52 0
199 31 23 33 0
200 52 40 50 0
201 43 52 56 0
202 33 23 26 0
203 50 40 39 0
204 56 52 66 0
205 50 39 51 0
206 66 52 63 0
207 56 66 75 0
208 51 39 41 0
209 63 52 50 0
210 41 39 31 0
211 51 41 54 0
212 41 31 33 0
213 54 41 44 0
214 41 33 44 0
215 44 33 36 0
216 65 80 76 0
217 65 76 62 0
218 42 55 49 0
219 42 49 37 0
220 35 43 46 0
221 46 43 56 0
222 46 56 59 0
223 59 56 75 0
224 46 59 53 0
225 53 59 67 0
226 67 59 71 0
227 71 59 75 0
228 67 71 81 0
229 81 71 89 0
230 67 81 72 0
231 72 81 91 0
232 91 81 100 0
233 100 81 89 0
234 91 100 112 0
235 112 100 120 0
236 120 100 114 0
237 112 120 132 0
238 114 100 89 0
239 120 114 133 0
240 132 120 141 0
241 133 114 125 0
242 120 133 141 0
243 125 114 110 0
244 125 110 137 0
245 133 125 146 0
246 125 137 146 0
247 144 136 122 0
248 144 122 132 0
249 132 122 112 0
250 140 165 160 0
251 66 63 82 0
252 236 240 241 0
253 241 240 244 0
254 134 154 137 0
255 134 137 115 0
256 50 51 60 0
257 60 51 64 0
258 50 60 63 0
259 64 51 54 0
260 60 64 74 0
261 63 60 73 0
262 74 64 79 0
263 60 74 73 0
264 79 64 68 0
265 74 79 93 0
266 73 74 90 0
267 68 64 54 0
268 93 79 95 0
269 73 90 88 0
270 95 79 83 0
271 93 95 111 0
272 88 90 106 0
273 83 79 68 0
274 111 95 116 0
275 106 90 108 0
276 88 106 101 0
277 116 95 102 0
278 111 116 130 0
279 108 90 93 0
280 102 95 83 0
281 116 102 121 0
282 130 116 135 0
283 93 90 74 0
284 102 83 87 0
285 121 102 109 0
286 135 116 121 0
287 87 83 70 0
288 102 87 109 0
289 70 83 68 0
290 87 70 77 0
291 109 87 96 0
292 70 68 57 0
293 57 68 54 0
294 70 57 61 0
295 57 54 44 0
296 57 44 47 0
297 108 93 111 0
298 108 111 127 0
299 127 111 130 0
300 108 127 123 0
301 127 130 148 0
302 123 127 145 0
303 148 130 152 0
304 127 148 145 0
305 152 130 135 0
306 145 148 167 0
307 152 135 159 0
308 167 148 174 0
309 145 167 162 0
310 159 135 142 0
311 152 159 173 0
312 174 148 152 0
313 162 167 187 0
314 142 135 121 0
315 173 159 196 0
316 187 167 192 0
317 162 187 181 0
318 142 121 129 0
319 192 167 174 0
320 187 192 211 0
321 181 187 203 0
322 129 121 109 0
323 192 174 205 0
324 203 187 211 0
325 181 203 204 0
326 204 203 223 0
327 181 204 184 0
328 184 204 208 0
329 181 184 163 0
330 208 204 219 0
331 184 208 188 0
332 163 184 168 0
333 181 163 162 0
334 219 204 223 0
335 184 188 168 0
336 163 168 147 0
337 162 163 143 0
338 147 168 150 0
339 163 147 143 0
340 143 147 124 0
341 162 143 145 0
342 124 147 126 0
343 143 124 123 0
344 123 124 106 0
345 106 124 101 0
346 145 143 123 0
347 126 147 150 0
348 123 106 108 0
349 126 150 134 0
350 152 173 174 0
351 174 198 205 0
352 63 73 82 0
353 73 88 82 0
354 129 109 118 0
355 124 126 101 0
356 240 239 243 0
357 209 188 208 0
358 209 208 221 0
359 221 208 219 0
360 209 221 222 0
361 222 221 230 0
362 209 222 214 0
363 230 221 235 0
364 222 230 232 0
365 232 230 237 0
366 222 232 225 0
367 142 129 151 0
368 151 129 139 0
369 142 151 166 0
370 166 151 186 0
371 142 166 159 0
372 159 166 196 0
373 94 107 97 0
374 226 231 234 0
375 226 234 229 0
376 128 140 119 0
377 193 201 212 0
378 193 212 206 0
379 126 134 115 0
380 115 137 110 0
381 221 219 235 0
382 244 240 243 0
383 97 78 76 0
384 36 47 44 0
385 223 228 219 0
386 5 3 7 0
387 198 174 173 0
388 37 27 29 0
389 96 118 109 0
390 242 243 239 0
391 183 160 185 0
392 6 10 8 0
393 190 196 166 0
394 190 166 186 0
395 62 48 49 0
396 61 77 70 0
397 235 237 230 0
398 1 2 3 0
399 211 217 203 0
400 18 11 16 0
401 139 161 151 0
402 245 244 243 0
403 138 119 140 0
404 17 26 23 0
405 78 62 76 0
406 47 61 57 0
407 228 235 219 0
408 205 211 192 0
409 118 139 129 0
410 160 138 140 0
411 10 17 14 0
412 196 198 173 0
413 48 37 49 0
414 77 96 87 0
415 237 242 232 0
416 2 6 4 0
417 217 223 203 0
418 11 5 7 0
419 119 97 107 0
420 26 36 33 0
421 89 110 114 0
422 101 82 88 0
423 82 75 66 0
424 115 101 126 0
425 75 89 71 0
426 165 185 160 0
427 151 161 186 0
428 2 4 3 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
2 1 0 
2 2 0 
6 1 0 
6 2 0 
10 1 0 
10 2 0 
17 1 0 
17 2 0 
26 1 0 
26 2 0 
36 1 0 
36 2 0 
47 1 0 
47 2 0 
61 1 0 
61 2 0 
77 1 0 
77 2 0 
96 1 0 
96 2 0 
118 1 0 
118 2 0 
139 1 0 
139 2 0 
161 1 0 
161 2 0 
185 1 0 
185 2 0 
189 1 0 
189 2 0 
197 1 0 
197 2 0 
199 1 0 
199 2 0 
206 1 0 
206 2 0 
212 1 0 
212 2 0 
216 1 0 
216 2 0 
224 1 0 
224 2 0 
229 1 0 
229 2 0 
234 1 0 
234 2 0 
238 1 0 
238 2 0 
241 1 0 
241 2 0 
244 1 0 
244 2 0 
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
1 1 1 
3 1 1 
5 1 1 
11 1 1 
18 1 1 
27 1 1 
37 1 1 
48 1 1 
62 1 1 
78 1 1 
97 1 1 
119 1 1 
138 1 1 
160 1 1 
182 1 -1 
183 1 1 
186 1 -1 
190 1 -1 
196 1 -1 
198 1 -1 
205 1 -1 
211 1 -1 
217 1 -1 
223 1 -1 
228 1 -1 
235 1 -1 
237 1 -1 
242 1 -1 
243 1 -1 
245 1 -1 
];

%% Volumetric Force
% Element        Dim                Force_Dim

Vol_force = [
];

%% Group Elements
% Element        Group_num

Group = [
];

%% Initial Holes
% Elements that are considered holes initially
% Element

Initial_holes = [
];

%% Boundary Elements
% Elements that can not be removed
% Element

Boundary_elements = [
];

%% Micro gauss post
%
% Element

Micro_gauss_post = [
];


%% Micro Slave-Master
% Nodes that are Slaves
% Nodes             Value (1-Slave,0-Master)

Micro_slave = [
];

%% Nodes solid
% Nodes that must remain 
% Nodes

nodesolid = unique(pointload_complete(:,1));

%% External border Elements
% Detect the elements that define the edge of the domain
% Element               Node(1)           Node(2)

External_border_elements = [
];

%% External border Nodes
% Detect the nodes that define the edge of the domain
% Node

External_border_nodes = [
];

%% Materials
% Materials that have been used
% Material_Num              Mat_density        Young_Modulus        Poisson

Materials = [
];
