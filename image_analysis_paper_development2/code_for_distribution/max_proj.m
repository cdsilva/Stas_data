function image1 = max_proj(zstack)

image1 = max(zstack, [], ndims(zstack));